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
h = cosmos.H0.to(u.km/u.s/u.Mpc).value / 100

### ========================================================== ###

# load the Planck kappa map
mask = hp.read_map('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/data/Planck/mask/mask.fits')
dat = hp.read_alm('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/data/Planck/MV/dat_klm.fits')
image = hp.sphtfunc.alm2map(dat, nside=2048, pol=False)
image_masked = hp.ma(image)
image_masked.mask = np.logical_not(mask)
print('finish loading Planck kappa map.')

### ========================================================== ###

quasar_cata = fits.getdata('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/data/SDSS16Q/DR16Q_v4.fits')
quasar_cata = quasar_cata[quasar_cata['Z'] > 1.5]
quasar_cata = quasar_cata[quasar_cata['Z'] < 2]


# TODO: read other file or do some cut
Nquas = len(quasar_cata)

c = coo.SkyCoord(ra=quasar_cata['RA']*u.degree, dec=quasar_cata['DEC']*u.degree)

l = c.galactic.l.to(u.rad).value
b = c.galactic.b.to(u.rad).value
pos = hp.ang2vec(theta=np.pi/2-b, phi=l)
# position of each quasar
w_l = np.ones(Nquas)
# weighting of each quasar
z_l = quasar_cata['Z']
# redshift of each quasar
print('finish loading quasar catalogue')

### ========================================================== ###
# random point sample

# data = np.loadtxt('./random_sample_theta_phi_z')
# pos = hp.ang2vec(theta=data[:,0], phi=data[:,1])
# z_l = data[:,2]
# Nquas = len(data)
# w_l = np.ones(Nquas)
# print('finish loading random sample')

### ========================================================== ###


Npro = 70
Njack = 100
Nbins = 15

r_bins = np.geomspace(3, 100, Nbins+1)        # unit: cMpc/h

sigma_frac = (const.c*const.c/(4*np.pi*const.G)).to(u.Msun*u.Mpc/u.pc/u.pc).value / h
z_s = 1100      # redshift of CMB
chi_s = coo.Distance(z=z_s, cosmology=cosmos).to(u.Mpc).value/(1+z_s)         # comoving distance of CMB

sender, receiver = Pipe(duplex=True)


def stack_single_sample(arg):
    global sender
    """
    To stack around a single position
    Estimator
    \Sigma_l(r) = (\sum_p w_lp \kappa_p \Sigma_c)/(\sum_p w_lp)
    w_lp = w_l * Sigma_c^(-2)
    \Sigma_c = c^c/(4\pi G) \chi_s/(\chi_l (\chi_s-\chi_l) (1+z_l)) 
    """
    vec = arg[0:3]
    w_l = arg[3]
    z_l = arg[4]
    values = np.zeros(Nbins)
    weights = np.zeros(Nbins)
    
    d_A = cosmos.angular_diameter_distance(z_l).to(u.Mpc).value
    chi_l = d_A*(1+z_l)

    theta_bins = r_bins/h/(1+z_l)/d_A
    # convert from cMpc to theta

    for i in range(Nbins):
        idx_in = hp.query_disc(nside=2048, vec=vec, radius=theta_bins[i])
        idx_out = hp.query_disc(nside=2048, vec=vec, radius=theta_bins[i+1])
        idx = np.zeros(len(image_masked), bool)
        idx[idx_out] = True
        idx[idx_in] = False
        # select pixels between bin[i] to bin[i+1]

        choose = image_masked[idx]
        kappa = choose[np.logical_not(choose.mask)].data
        
        Sigma_c = sigma_frac*chi_s/(chi_l*(chi_s-chi_l)*(1+z_l))*np.ones_like(kappa)        # unit: M_sun/pc^2 h
        # Sigma_c = np.ones_like(kappa)

        w_lp = w_l/Sigma_c/Sigma_c

        weights[i] = np.sum(w_lp)
        if len(w_lp) == 0:
            values[i] = np.nan
        else:
            values[i] = np.sum(w_lp*Sigma_c*kappa)/weights[i]
    if not sender is None:
        sender.send(1)

    return values, weights

# def boostrap_err(value, weight, Nboostrap):
#     boostrap_value = np.zeros(Nboostrap)
#     for i in range(Nboostrap):
#         idx = np.random.randint(0, len(value), len(value))
#         _value = value[idx]
#         _weight = weight[idx]
#         if len(_weight) == 0:
#             boostrap_value[i] = np.nan
#         else:
#             boostrap_value[i] = np.nansum(_value*_weight)/_weight.sum()
#     return np.nanstd(boostrap_value)

def jackknife(value, weight, Njack):
    # should not have nan.
    N = len(value)
    Nsub = int(N/Njack)
    val_jackknife = np.zeros(Njack)
    for i in range(Njack):
        idx = int(i*N/Njack)
        # idx~idx+Nsub will be subtracted
        index = np.ones(N, bool)
        index[idx:idx+Nsub] = False
        val = value[index]
        wei = weight[index]
        val_jackknife[i] = (val*wei).sum()/wei.sum()
    return val_jackknife.std()*np.sqrt(Njack)

def progress_bar(receiver, rate=10):
    # rate is the update rate of the progress bar.
    rate = max(1, int(rate))
    with tqdm.tqdm(total=Nquas) as pbar:
        count = 0
        while True:
            if receiver.recv():
                count += 1
                if count % rate == 0:
                    pbar.update(rate)
            else:
                break
        pbar.close()

mypool = Pool(processes=Npro)

args = np.hstack((pos, w_l.reshape((-1, 1)), z_l.reshape(-1, 1)))
assert args.shape == (Nquas, 5)
print('start calculation...')

progress = Process(target=progress_bar, args=(receiver,100))
progress.start()

res_async = mypool.map_async(stack_single_sample, args, chunksize=200)

result = res_async.get()
sender.send(0)

progress.join()

values = np.vstack(tuple(item[0] for item in result))
weight = np.vstack(tuple(item[1] for item in result))
# values, weight should have shape (Nquas, Nbins)
assert values.shape==weight.shape==(Nquas, Nbins)

print('finish calculating')
print('completeness: {}'.format(1-np.sum(np.isnan(values).flat)/(Nquas*Nbins)))

np.save('./calculation_data/result_r=3_100_15.npy', (values, weight))

print('start stacking and jackknife...')
sigma = np.zeros(Nbins)
sigma_err = np.zeros(Nbins)
for i in range(Nbins):
    w = weight[:,i]
    v = values[:,i]
    idx = np.logical_not(np.isnan(v))

    v = v[idx]
    w = w[idx]
    w /= w.mean()
    sigma[i] = np.sum(w*v)/w.sum()
    sigma_err[i] = jackknife(v, w, Njack)

import pandas as pd
result = pd.DataFrame({'r': (r_bins[1:]+r_bins[:-1])/2, 'sigma': sigma, 'sigma_err': sigma_err})
result.to_csv('./results/result_r=3_100_15')
