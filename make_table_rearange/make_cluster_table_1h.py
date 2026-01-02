import numpy as np
import pickle
from Corrfunc.theory import DD, DDrppi
from Corrfunc.utils import convert_3d_counts_to_cf, convert_rp_pi_counts_to_wp
import tqdm
import matplotlib.pyplot as plt
with open('./displacement.bin', 'rb') as f:
    displacement = pickle.load(f)
nthread = 40

sep_min = 0.01
sep_max = 100
sep_N = 50

pimax = 100                         # this should match the result from auto corr
boxsize = 2500

rp_bins = np.geomspace(sep_min, sep_max, sep_N+1)
rp = (rp_bins[1:]*rp_bins[:-1])**0.5

## 1h s-s term
NR = 1000_000_000
n = NR / (boxsize**3)

wp_array = np.zeros((len(displacement), sep_N))

counts = []
for i in tqdm.trange(len(displacement)):
    sat = displacement[i]
    ND = len(sat)
    D1D2 = DDrppi(1, nthread, pimax, rp_bins, X1=sat[:,0].T + boxsize/2, Y1=sat[:,1].T + boxsize/2, Z1=sat[:,2].T + boxsize/2, periodic=False, boxsize=boxsize, verbose=False)

    counts.append(D1D2)

import pickle
with open('bin/1h_ss_count.bin', 'wb') as f:
    pickle.dump(counts, f)