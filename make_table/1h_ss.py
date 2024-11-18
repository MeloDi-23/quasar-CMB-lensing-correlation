import numpy as np
import pickle
from Corrfunc.theory import DD, DDrppi
from Corrfunc.utils import convert_3d_counts_to_cf, convert_rp_pi_counts_to_wp
import tqdm
# 1h s-s term
with open('./displacement.bin', 'rb') as f:
    displacement = pickle.load(f)
nthread = 40

sep_min = 0.1
sep_max = 100
sep_N = 30

pimax = 100                         # this should match the result from auto corr
boxsize = 2500

rp_bins = np.geomspace(sep_min, sep_max, sep_N+1)
NR = 1000_000_000
n = NR / (boxsize**3)

wp_array = np.zeros((len(displacement), sep_N))

counts = []
for i in tqdm.trange(len(displacement)):
    sat = displacement[i]
    ND = len(sat)
    D1D2 = DDrppi(1, nthread, pimax, rp_bins, X1=sat[:,0].T + boxsize/2, Y1=sat[:,1].T + boxsize/2, Z1=sat[:,2].T + boxsize/2, periodic=False, boxsize=boxsize, verbose=False)

    vol = np.pi*(D1D2['rmax']**2 - D1D2['rmin']**2)*1           # delta r pi = 1

    counts.append(D1D2)


    RR = D1D2.copy()

    RR['npairs'] = (NR*vol*n*2).astype('uint64')            # this has been cross checked to be correct
    
    D1R2 = RR.copy()
    D1R2['npairs'] = ND*vol*n*2

    wp_array[i] = convert_rp_pi_counts_to_wp(ND, ND, NR, NR, D1D2, D1R2, D1R2, RR, pimax=pimax, nrpbins=sep_N)

import pickle
with open('1h_ss_count.bin', 'wb') as f:
    pickle.dump(counts, f)
np.save('1h_ss_wp', wp_array)