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
Nside = 1024
mask = hp.read_map('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/data/Planck/mask/mask.fits')
dat = hp.read_alm('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/data/Planck/MV/dat_klm.fits')
image = hp.sphtfunc.alm2map(dat, nside=Nside, pol=False)
mask = hp.ud_grade(mask, Nside)
image_masked = hp.ma(image)
image_masked.mask = np.logical_not(mask>0.5)
print('finish loading Planck kappa map.')

### ========================================================== ###

cmass = np.load('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/cmass_z_cut.npy')
Nquas = len(cmass)
c = coo.SkyCoord(ra=cmass['ra']*u.degree, dec=cmass['dec']*u.degree)

l = c.galactic.l.to(u.rad).value
b = c.galactic.b.to(u.rad).value
pos = hp.ang2vec(theta=np.pi/2-b, phi=l)
# position of each quasar
w_l = cmass['w']
# weighting of each quasar
z_l = cmass['z']
# redshift of each quasar
print('finish loading cmass catalogue')

### ========================================================== ###
# random point sample
# cmass = np.load('/uufs/astro.utah.edu/common/home/u6060319/quasar-CMBlening/cmass_random.npy')

# Nquas = len(cmass)

# c = coo.SkyCoord(ra=cmass['ra']*u.degree, dec=cmass['dec']*u.degree)

# l = c.galactic.l.to(u.rad).value
# b = c.galactic.b.to(u.rad).value
# pos = hp.ang2vec(theta=np.pi/2-b, phi=l)
# # position of each quasar
# w_l = cmass['w']
# # weighting of each quasar
# z_l = cmass['z']
# # redshift of each quasar

# print('finish loading random sample')

### ========================================================== ###

Npro = 60
Njack = 100
Nbins = 15
sep_min = 0.5
sep_max = 100
name = 'cmass'
r_bins = np.geomspace(sep_min, sep_max, Nbins+1)        # unit: cMpc/h

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
        idx_in = hp.query_disc(nside=Nside, vec=vec, radius=theta_bins[i])
        idx_out = hp.query_disc(nside=Nside, vec=vec, radius=theta_bins[i+1])
        idx = np.zeros(len(image_masked), bool)
        idx[idx_out] = True
        idx[idx_in] = False
        # select pixels between bin[i] to bin[i+1]

        choose = image_masked[idx]
        kappa = choose[np.logical_not(choose.mask)].data

        Sigma_c = sigma_frac*chi_s/(chi_l*(chi_s-chi_l)*(1+z_l))*np.ones_like(kappa)        # unit: M_sun/pc^2 h

        w_lp = w_l/Sigma_c/Sigma_c

        weights[i] = np.sum(w_lp)
        if len(w_lp) == 0:
            values[i] = np.nan
        else:
            values[i] = np.sum(w_lp*Sigma_c*kappa)/weights[i]
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

args = np.hstack((pos, w_l.reshape((-1, 1)), z_l.reshape(-1, 1)))
assert args.shape == (Nquas, 5)
print('start calculation...')

progress = Process(target=progress_bar, args=(receiver,100))
progress.start()

result = mypool.map(stack_single_sample, args, chunksize=200)
# result = res_async.get()

sender.send(0)
progress.join()

values = np.vstack(tuple(item[0] for item in result))
weight = np.vstack(tuple(item[1] for item in result))
assert values.shape==weight.shape==(Nquas, Nbins)

print('finish calculating')
print('completeness: {}'.format(1-np.sum(np.isnan(values).flat)/(Nquas*Nbins)))

np.save(f'./calculation_data/result_r={sep_min}_{sep_max}_{Nbins}_{name}_low.npy', (values, weight))

# print('start stacking and jackknife...')
# sigma = np.zeros(Nbins)
# sigma_err = np.zeros(Nbins)

# for i in range(Nbins):
#     w = weight[:,i]
#     v = values[:,i]
#     idx = np.logical_not(np.isnan(v))

#     v = v[idx]
#     w = w[idx]
#     w /= w.mean()
#     sigma[i] = np.sum(w*v)/w.sum()


# def weight_nan_mean(value, weight):
#     mean_val = np.zeros(value.shape[1])
#     for i in range(len(mean_val)):
#         idx = np.logical_not(np.isnan(value[:,i]))
#         val = value[:,i][idx]
#         wei = weight[:,i][idx]
#         mean_val[i] = (val*wei).sum()/wei.sum()
#     return mean_val

# Njack = 100
# Nsub = int(Nquas/Njack)
# shuffled_index = np.arange(Nquas)
# np.random.shuffle(shuffled_index)

# pix = hp.vec2pix(3, pos[:,0], pos[:,1], pos[:,2])
# def jackknife_resample(i):
#     idx = pix != i
#     return weight_nan_mean(values[idx], weight[idx])

# mypool = Pool(Npro)
# result = mypool.map(jackknife_resample, np.unique(pix))
# jackknife_value = np.vstack(result)

# # jackknife_nside = hp.pixelfunc.get_min_valid_nside(Njack)
# # jackknife_npix = hp.nside2npix(jackknife_nside)
# # pix = hp.vec2pix(jackknife_nside, pos[:,0], pos[:,1], pos[:,2])
# # jackknife_value = np.zeros((Nbins, Njack))
# # for i in range(Nbins):
# #     for j in range(Njack):
# #         include = pix != np.random.randint(jackknife_npix)
# #         notnan = np.logical_not(np.isnan(v))
# #         idx = np.logical_and(include, notnan)
# #         v = v[idx]
# #         w = w[idx]
# #         w /= w.mean()
# #         jackknife_value[i,j] = np.sum(w*v)/w.sum()

# sigma_err = jackknife_value.std(axis=0)*np.sqrt(Njack)

# import pandas as pd
# result = pd.DataFrame({'r': (r_bins[1:]+r_bins[:-1])/2, 'sigma': sigma, 'sigma_err': sigma_err})
# result.to_csv(f'./results/result_r={sep_min}_{sep_max}_{Nbins}_{name}_low')
