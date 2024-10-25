import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import Planck18 as cosmos
import astropy.coordinates as coo
from Corrfunc.mocks import DDrppi_mocks, DDtheta_mocks
from Corrfunc.io import read_catalog
from Corrfunc.theory import DDrppi, DD
from Corrfunc.utils import convert_rp_pi_counts_to_wp, convert_3d_counts_to_cf
from tqdm import tqdm
import multiprocessing as mp
from itertools import repeat
import healpy as hp

h = cosmos.H0.value/100

cmass = np.load('../catalogue/cmass_z_cut.npy')
random = np.load('../catalogue/random_cmass.npy')
Nbins = 16
rp_bin = np.geomspace(0.1, 150, Nbins+1)

Nd = len(cmass)
Nr = len(random)

pimax = 160

dis_cov_g = cosmos.angular_diameter_distance(cmass['z']).to(u.Mpc).value*(1+cmass['z'])*h
dis_cov_r = cosmos.angular_diameter_distance(random['z']).to(u.Mpc).value*(1+random['z'])*h

Nside_jack = 6
c = coo.SkyCoord(ra=cmass['ra']*u.degree, dec=cmass['dec']*u.degree)

l = c.galactic.l.to(u.rad).value
b = c.galactic.b.to(u.rad).value
npix = hp.ang2pix(Nside_jack, theta=np.pi/2-b, phi=l)
c = coo.SkyCoord(ra=random['ra']*u.degree, dec=random['dec']*u.degree)

l = c.galactic.l.to(u.rad).value
b = c.galactic.b.to(u.rad).value
npix_r = hp.ang2pix(Nside_jack, theta=np.pi/2-b, phi=l)

pix = np.unique(npix)

def resample(p, sender):
    sub_sample = npix != p
    sub_sample_r = npix_r != p

    cmass_sub = cmass[sub_sample]
    random_sub = random[sub_sample_r]
    dis_sub_g = dis_cov_g[sub_sample]
    dis_sub_r = dis_cov_r[sub_sample_r]

    dd = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                   RA1=cmass_sub['ra'], DEC1=cmass_sub['dec'], CZ1=dis_sub_g, weights1=cmass_sub['w'], is_comoving_dist=True, weight_type='pair_product')
    rr = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin,
                    RA1=random_sub['ra'], DEC1=random_sub['dec'], CZ1=dis_sub_r, weights1=random_sub['w'], is_comoving_dist=True, weight_type='pair_product')
    dr = DDrppi_mocks(
        autocorr=False, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin, 
        RA1=cmass_sub['ra'], DEC1=cmass_sub['dec'], CZ1=dis_sub_g, weights1=cmass_sub['w'], 
        RA2=random_sub['ra'], DEC2=random_sub['dec'], CZ2=dis_sub_r, weights2=random_sub['w'], 
        is_comoving_dist=True, weight_type='pair_product')

    Nd = len(cmass_sub)
    Nr = len(random_sub)
    sender.send(1)
    return convert_rp_pi_counts_to_wp(Nd, Nd, Nr, Nr, dd, dr, dr, rr, pimax=pimax, nrpbins=Nbins)

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
sender, recv = mp.Pipe()

pool = mp.Pool(20)
print('start')
bar = mp.Process(target=progress_bar, args=(len(pix), recv, 1))
bar.start()
w_arr = np.vstack(pool.starmap(resample, zip(pix, repeat(sender))))
sender.send(0)
bar.join()
print('end')

np.save('auto_corr_jackknife_cmass', w_arr)

