import numpy as np
from Corrfunc.theory import DD, DDrppi
from Corrfunc.utils import convert_3d_counts_to_cf, convert_rp_pi_counts_to_wp

import pickle


halos = np.load('main_halo.npy')
particle = np.loadtxt('./dm_xyz_snap_010_thin.dat')


logMh_m = 11.5
logMh_M = 14.72
logMh_bin = 0.02
logMh_N = int((logMh_M - logMh_m) / logMh_bin)

bins = np.arange(logMh_N+1) * logMh_bin + logMh_m
res = np.digitize(np.log10(halos['M_h']), bins)

sep_min = 0.01
sep_max = 100
sep_N = 50

boxsize = 2500

r_pbins = np.geomspace(sep_min, sep_max, sep_N+1)
nthread = 40
# res = 0, 1, ... N+1

result = []
import tqdm

for i in tqdm.tqdm(range(1, logMh_N+1)):
    halo1 = halos[res == i]
    if len(halo1) == 0:
        result.append(None)
        continue

    ret = DD(0, nthread, r_pbins, X1=halo1['X'], Y1=halo1['Y'], Z1=halo1['Z'], 
                X2=particle[:,0].T, Y2=particle[:,1].T, Z2=particle[:,2].T, periodic=True, boxsize=boxsize, verbose=False)
    result.append(ret)


with open(f'satellite_distribution.bin', 'wb') as f:
    pickle.dump(result, f)
