import pycorr
# import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18 as cosmos
from astropy import units as u
from astropy.io import fits
from sklearn.neighbors import KDTree
from data_process import jackknife_label
h = cosmos.H0.value/100
print('program start', flush=True)
quasar = np.load('../catalogue/quasar_lss_all.npy')
random = np.load('../catalogue/random_quasar_lss_all.npy')
npix = jackknife_label(quasar, 10)
npix_r = jackknife_label(random, 10)
pix = np.unique(npix)

label = np.zeros_like(npix)
label_r = np.zeros_like(npix_r)
all_pix = np.unique(npix)
for i, p in enumerate(all_pix):
    label[npix==p] = i
    label_r[npix_r==p] = i
dis_cov_q = cosmos.comoving_distance(quasar['z']).to(u.Mpc).value*h
dis_cov_r = cosmos.comoving_distance(random['z']).to(u.Mpc).value*h
import astropy.units as u
from astropy.cosmology import Planck18 as cosmos
import astropy.coordinates as coo
h = cosmos.H0.value/100

Nbins = 15
rp_bin = np.geomspace(3, 100, Nbins+1)
pimax = 100
pos_q = np.vstack([quasar['ra'], quasar['dec'], dis_cov_q])
pos_r = np.vstack([random['ra'], random['dec'], dis_cov_r])
print('start calculation', flush=True)

correlation_func = pycorr.correlation_function.TwoPointCorrelationFunction(
    'rppi', (rp_bin, np.linspace(-pimax, pimax, 2*pimax+1, endpoint=True)),
    pos_q, pos_q, pos_r, pos_r, 
    data_weights1=quasar['w'],
    data_weights2=quasar['w'],
    data_samples1=label,
    data_samples2=label,
    randoms_samples1=label_r,
    randoms_samples2=label_r,
    randoms_weights1=random['w'],
    randoms_weights2=random['w'],
    estimator='landyszalay',
    position_type='rdd'
)

print('saving calculation', flush=True)

import pickle
with open('pycorr_1h.bin', 'wb') as f:
    pickle.dump(correlation_func, f)