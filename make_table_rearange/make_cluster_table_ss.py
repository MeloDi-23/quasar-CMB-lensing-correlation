import numpy as np
from Corrfunc.theory import DD, DDrppi
import pickle


def pair_count(halo1, halo2, autocorr, nthread, rp_bins, pimax, boxsize):
    if autocorr:
        D1D2_rppi = DDrppi(autocorr, nthread, pimax, rp_bins, X1=halo1['X'], Y1=halo1['Y'], Z1=halo1['Z'], periodic=True, boxsize=boxsize, verbose=False)
    else:
        D1D2_rppi = DDrppi(autocorr, nthread, pimax, rp_bins, X1=halo1['X'], Y1=halo1['Y'], Z1=halo1['Z'], X2=halo2['X'], Y2=halo2['Y'], Z2=halo2['Z'], periodic=True, boxsize=boxsize, verbose=False)
    return D1D2_rppi

halos = np.load('satellite_halo_zspace.npy')
logMh_m = 11.5
logMh_M = 14.72
logMh_bin = 0.02
logMh_N = int((logMh_M - logMh_m) / logMh_bin)

bins = np.arange(logMh_N+1) * logMh_bin + logMh_m
res = np.digitize(np.log10(halos['M_h']), bins)

sep_min = 0.01
sep_max = 100
sep_N = 50

pimax = 100                         # this should match the result from auto corr
boxsize = 2500

r_pbins = np.geomspace(sep_min, sep_max, sep_N+1)
nthread = 40
# res = 0, 1, ... N+1

result = {}
import tqdm

for i in tqdm.tqdm(range(1, logMh_N+1)):
    halo1 = halos[res == i]
    if len(halo1) == 0:
        continue
    for j in range(1, i):
        halo2 = halos[res == j]
        if len(halo2) == 0:
            continue
        ret = pair_count(halo1, halo2, 0, nthread, r_pbins, pimax, boxsize)
        result[(i,j)] = ret
    ret = pair_count(halo1, halo1, 1, nthread, r_pbins, pimax, boxsize)
    result[(i,i)] = ret

with open('bin/s_s_pair_count.bin', 'wb') as f:
    pickle.dump(result, f)
