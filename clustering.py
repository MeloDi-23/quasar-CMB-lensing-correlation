import numpy as np
import matplotlib.pyplot as plt
from data_process import jackknife_label
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
cosmos = FlatLambdaCDM(H0=67.77, Om0=0.307115, Ob0=0.048206, Tcmb0=2.7255)
import astropy.coordinates as coo
h = cosmos.H0.value/100

from Corrfunc.mocks import DDrppi_mocks
from my_util import convert_rp_pi_counts_to_wp
import multiprocessing as mp
from tqdm import tqdm

import os
import sys


if len(sys.argv) < 3:
    sys.exit(-1)
file_data = sys.argv[1]
file_random = sys.argv[2]
data_name = os.path.splitext(os.path.split(file_data)[1])[0]

print(data_name)
# sys.exit(0)

file_save = os.path.join('./auto_correlation', data_name+'_auto.npy')
print(f'save to {file_save}')

Nbins = 15
rp_min = 3
rp_max = 100

rp_bin_cross = np.geomspace(rp_min, rp_max, Nbins+1)
r_p_cross = (rp_bin_cross[:-1]*rp_bin_cross[1:])**0.5

sep = (np.log10(rp_max) - np.log10(rp_min)) / Nbins
rp_bin_auto = 10**(np.arange(-10, Nbins+1)*sep + np.log10(rp_min))
# rp_bin_auto = np.geomspace(rp_min, rp_max, Nbins+1)
r_p_auto = (rp_bin_auto[:-1]*rp_bin_auto[1:])**0.5

Nbins = len(rp_bin_auto) - 1

Nside_jack = 10
quasar_cata = np.load(file_data)
random_sample = np.load(file_random)
pix = jackknife_label(quasar_cata, Nside_jack)
pix_r = jackknife_label(random_sample, Nside_jack)

Nd = len(quasar_cata)
Nr = len(random_sample)
dis_cov_q = cosmos.angular_diameter_distance(quasar_cata['z']).to(u.Mpc).value*(1+quasar_cata['z'])*h
dis_cov_r = cosmos.angular_diameter_distance(random_sample['z']).to(u.Mpc).value*(1+random_sample['z'])*h

pimax = 100

npix = pix
npix_r =  pix_r

pix = np.unique(npix)
def resample(p, sender):
    sub_sample = npix != p
    sub_sample_r = npix_r != p

    quasar_sub = quasar_cata[sub_sample]
    random_sub = random_sample[sub_sample_r]

    dd = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin_auto,
                   RA1=quasar_sub['ra'], DEC1=quasar_sub['dec'], CZ1=dis_cov_q[sub_sample], weights1=quasar_sub['w'], is_comoving_dist=True, weight_type='pair_product')
    rr = DDrppi_mocks(autocorr=True, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin_auto,
                    RA1=random_sub['ra'], DEC1=random_sub['dec'], CZ1=dis_cov_r[sub_sample_r], weights1=random_sub['w'], is_comoving_dist=True, weight_type='pair_product')
    dr = DDrppi_mocks(
        autocorr=False, cosmology=1, nthreads=50, pimax=pimax, binfile=rp_bin_auto, 
        RA1=quasar_sub['ra'], DEC1=quasar_sub['dec'], CZ1=dis_cov_q[sub_sample], weights1=quasar_sub['w'], 
        RA2=random_sub['ra'], DEC2=random_sub['dec'], CZ2=dis_cov_r[sub_sample_r], weights2=random_sub['w'], 
        is_comoving_dist=True, weight_type='pair_product')

    Nd = quasar_sub['w'].sum()
    Nr = random_sub['w'].sum()
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

from itertools import repeat

sender, recv = mp.Pipe()

pool = mp.Pool(5)

bar = mp.Process(target=progress_bar, args=(len(pix), recv, 1))
bar.start()
jack_auto = np.vstack(list(pool.starmap(resample, zip(pix, repeat(sender)))))
sender.send(0)
bar.join()

np.save(file_save, jack_auto)