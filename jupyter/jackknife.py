"""
Used to calculate jackknife of things you want.
copy from .ipynb file.
"""

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from tqdm import tqdm
import astropy.units as u
from astropy.cosmology import Planck18 as cosmos
import astropy.coordinates as coo
h = cosmos.H0.value/100
from Corrfunc.mocks import DDrppi_mocks, DDtheta_mocks
from Corrfunc.io import read_catalog
from Corrfunc.theory import DDrppi, DD
from my_util import convert_rp_pi_counts_to_wp, convert_3d_counts_to_cf
from astropy.io import fits
from sklearn.neighbors import KDTree
quasar = np.load('../catalogue/quasar_lss_all.npy')
random = np.load('../catalogue/random_quasar_lss_all.npy')

Nd = len(quasar)
Nr = len(random)
dis_cov_q = cosmos.angular_diameter_distance(quasar['z']).to(u.Mpc).value*(1+quasar['z'])*h
dis_cov_r = cosmos.angular_diameter_distance(random['z']).to(u.Mpc).value*(1+random['z'])*h
quasar_SDSS = fits.getdata('/uufs/chpc.utah.edu/common/home/astro/zheng/hd/data/SDSS16Q/DR16Q_v4.fits')
tree = KDTree(np.c_[quasar_SDSS['RA'], quasar_SDSS['DEC']], metric='euclidean')
que = tree.query(np.c_[quasar['ra'], quasar['dec']])
valid = que[0].flatten() < 5/3600
index = que[1].flatten()

M_I = quasar_SDSS['M_I'][index]
M_I[~valid] = np.nan
bins = np.linspace(0.8, 2.2, 30)            # the z cut applied to quasar lss all
result = np.digitize(quasar['z'], bins)
kind = np.zeros(len(quasar), int)
middles = []
high = []
low = []
for i in range(1, 30):
    index = np.where(result == i)[0]
    M = M_I[index]
    middle = np.percentile(M[~np.isnan(M)], 50)
    high.append(index[M <= middle])
    low.append(index[M >= middle])

index_h = np.concatenate(high)
index_l = np.concatenate(low)
Nbins = 15
rp_min = 3
rp_max = 100
rp_bin = np.geomspace(rp_min, rp_max, Nbins+1)
r_p = (rp_bin[:-1]*rp_bin[1:])**0.5

z = 1.69
h = cosmos.H0.to(u.km/u.s/u.Mpc).value / 100
d_A = cosmos.angular_diameter_distance(z).to(u.Mpc).value
chi_l = d_A*(1+z)

theta_bins = rp_bin/h/(1+z)/d_A
theta_bins_deg = np.rad2deg(theta_bins)

pimax = 100


with open('label.bin', 'rb') as f:
    npix = np.load(f)
    npix_r = np.load(f)
pix = np.unique(npix)

_npix = npix
_npix_r = npix_r
_quasar = quasar
_random = random
_dis_cov_q = dis_cov_q
_dis_cov_r = dis_cov_r


def progress_bar(Ntotal, receiver, rate=10):
    # rate is the update rate of the progress bar.
    rate = max(1, int(rate))
    with tqdm(total=int(Ntotal/rate)) as pbar:
        count = 0
        while True:
            if receiver.recv():
                count += 1
                if count % rate == 0:
                    pbar.update(1)
            else:
                break
        pbar.close()
from itertools import repeat
sender, recv = mp.Pipe()



for i, index in zip(['high', 'low'], [index_h, index_l]):
    _npix = npix[index]
    _quasar = quasar[index]
    _dis_cov_q = dis_cov_q[index]
    bar = mp.Process(target=progress_bar, args=(len(pix), recv, 1))
    bar.start()
    def resample(p, sender):
        sub_sample = _npix != p
        sub_sample_r = _npix_r != p

        quasar_sub = _quasar[sub_sample]
        random_sub = _random[sub_sample_r]

        dd = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                    RA1=quasar_sub['ra'], DEC1=quasar_sub['dec'], CZ1=_dis_cov_q[sub_sample], weights1=quasar_sub['w'], 
                    is_comoving_dist=True, weight_type='pair_product')
        rr = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                    RA1=random_sub['ra'], DEC1=random_sub['dec'], CZ1=_dis_cov_r[sub_sample_r], weights1=random_sub['w'], 
                    is_comoving_dist=True, weight_type='pair_product')
        dr = DDrppi_mocks(
            autocorr=False, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin, 
            RA1=quasar_sub['ra'], DEC1=quasar_sub['dec'], CZ1=_dis_cov_q[sub_sample], weights1=quasar_sub['w'], 
            RA2=random_sub['ra'], DEC2=random_sub['dec'], CZ2=_dis_cov_r[sub_sample_r], weights2=random_sub['w'], 
            is_comoving_dist=True, weight_type='pair_product')

        Nd = quasar_sub['w'].sum()
        Nr = len(random_sub)
        sender.send(1)
        return convert_rp_pi_counts_to_wp(Nd, Nd, Nr, Nr, dd, dr, dr, rr, pimax=pimax, nrpbins=Nbins)
    pool = mp.Pool(5)
    w_arr = np.vstack(list(pool.starmap(resample, zip(pix, repeat(sender)))))
    sender.send(0)
    bar.join()
    np.save(f'auto_corr_jackknife_kmeans_100_{i}', w_arr)