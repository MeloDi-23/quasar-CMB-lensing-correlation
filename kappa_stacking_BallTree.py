import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
from multiprocessing import Pool, Pipe, Process
import tqdm
from astropy.cosmology import Planck18 as cosmos
import astropy.coordinates as coo
# from sklearn.neighbors import BallTree
h = cosmos.H0.to(u.km/u.s/u.Mpc).value / 100
import argparse
"""
Use BallTree to calculate the result, prove to be consistent with hp.query, but seems faster.
Supports command line and parameters.
Must use the galaxy catalogue of format: .npy, with ra, dec, z, w(weight)
Must use the kappa map of format: .npy, with l, b, kappa
"""
### ========================================================== ###
parser = argparse.ArgumentParser('kappa stacking (balltree)')
parser.add_argument('file', help='input catalogue file')
parser.add_argument('CMB_file', help='input CMB kappa map file')
parser.add_argument('--sep_min', '-m', default=3, help='min seperation [cMpc/h]')
parser.add_argument('--sep_max', '-M', default=100, help='max seperation [cMpc/h]')
parser.add_argument('--nbins', '-N', default=10, type=int, help='number of bins')
parser.add_argument('--name', '-n', help='output name')
parser.add_argument('--frame', '-f', help='frame of coordinats for catalogue file', default='icrs', choices=['galactic', 'icrs'])
parser.add_argument('--unit', '-u', help='unit of catalogue file', default='deg')
arg = parser.parse_args()

if not arg.name:
    import os
    file_name = os.path.split(arg.file)[1]
    arg.name = os.path.splitext(file_name)[0] + '_' + os.path.splitext(os.path.split(arg.CMB_file)[1])[0]
### ========================================================== ###

# load the Planck kappa map
print(f'loading {arg.CMB_file}...')
l_k, b_k, kappa = np.load(arg.CMB_file)

from sklearn.neighbors import BallTree
tree = BallTree(data=np.vstack((b_k, l_k)).T, leaf_size=5, metric='haversine')          # latitude + logtitude

print('finish loading Planck kappa map.')

### ========================================================== ###
print(f'loading {arg.file}...')
catalogue = np.load(arg.file)
Nquas = len(catalogue)
c = coo.SkyCoord(l=catalogue['l'], b=catalogue['b'], frame=arg.frame, unit=arg.unit)

l = c.galactic.l.to(u.rad).value
b = c.galactic.b.to(u.rad).value
# position of each quasar
w_l = catalogue['w']
# weighting of each quasar
z_l = catalogue['z']
# redshift of each quasar
print('finish loading catalogue.')

### ========================================================== ###

Npro = 60
Njack = 100
Nbins = arg.nbins
sep_min = float(arg.sep_min)
sep_max = float(arg.sep_max)
name = arg.name
r_bins = np.geomspace(sep_min, sep_max, Nbins+1)        # unit: cMpc/h

output = f'./calculation_data/result_r={arg.sep_min}_{arg.sep_max}_{Nbins}_{name}_tree.npy'
print(f'output path {output}')

sigma_frac = (const.c*const.c/(4*np.pi*const.G)).to(u.Msun*u.Mpc/u.pc/u.pc).value / h
z_s = 1100      # redshift of CMB
chi_s = coo.Distance(z=z_s, cosmology=cosmos).to(u.Mpc).value/(1+z_s)         # comoving distance of CMB

sender, receiver = Pipe()

def stack_single_sample(arg):
    global sender
    r"""
    To stack around a single position
    Estimator
    \Sigma_l(r) = (\sum_p w_lp \kappa_p \Sigma_c)/(\sum_p w_lp)
    w_lp = w_l * Sigma_c^(-2)
    \Sigma_c = c^c/(4\pi G) \chi_s/(\chi_l (\chi_s-\chi_l) (1+z_l))
    """
    ra, dec = arg[0], arg[1]
    w_l = arg[2]
    z_l = arg[3]
    values = np.zeros(Nbins)
    weights = np.zeros(Nbins)

    d_A = cosmos.angular_diameter_distance(z_l).to(u.Mpc).value
    chi_l = d_A*(1+z_l)

    theta_bins = r_bins/h/(1+z_l)/d_A
    # convert from cMpc to theta

    for i in range(Nbins):
        idx_in, idx_out = tree.query_radius([[dec, ra], [dec, ra]], [theta_bins[i], theta_bins[i+1]])
        idx = np.zeros(len(kappa), bool)
        idx[idx_out] = True
        idx[idx_in] = False
        # select pixels between bin[i] to bin[i+1]

        k = kappa[idx]

        Sigma_c = sigma_frac*chi_s/(chi_l*(chi_s-chi_l)*(1+z_l))*np.ones_like(k)        # unit: M_sun/pc^2 h

        w_lp = w_l/Sigma_c/Sigma_c

        weights[i] = np.sum(w_lp)
        if len(w_lp) == 0:
            values[i] = np.nan
        else:
            values[i] = np.sum(w_lp*Sigma_c*k)/weights[i]
    if not sender is None:
        sender.send(1)

    return values, weights

def progress_bar(receiver, rate=10):
    # rate is the update rate of the progress bar.
    rate = max(1, int(rate))
    with tqdm.tqdm(total=int(Nquas/rate)) as pbar:
        count = 0
        while True:
            if receiver.recv():
                count += 1
                if count % rate == 0:
                    pbar.update(1)
            else:
                break
        pbar.close()

mypool = Pool(processes=Npro)

args = np.vstack((l, b, w_l, z_l)).T
assert args.shape == (Nquas, 4)
print('start calculation...')

progress = Process(target=progress_bar, args=(receiver,100))
progress.start()

result = mypool.map(stack_single_sample, args, chunksize=200)

sender.send(0)
progress.join()

values = np.vstack(tuple(item[0] for item in result))
weight = np.vstack(tuple(item[1] for item in result))
assert values.shape==weight.shape==(Nquas, Nbins)

print('finish calculating.')
print('completeness: {}.'.format(1-np.sum(np.isnan(values))/(Nquas*Nbins)))

print(f'save to {output}...')
np.save(output, (values, weight))

