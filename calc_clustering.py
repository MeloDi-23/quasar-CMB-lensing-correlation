import numpy as np
from NN_correlation import two_point_counter_2d
from astropy.cosmology import Planck18 as cosmos
h = cosmos.H0.value/100

quasar = np.load('./catalogue/quasar_lss.npy')
random = np.load('./catalogue/random_quasar_lss.npy')

sep_min = 3
sep_max = 100
Nbins = 15
rp_bins = np.geomspace(sep_min, sep_max, Nbins+1)
rp = np.sqrt(rp_bins[1:]*rp_bins[:-1])
pimax = 40
Npibins = 40
pi_bins = np.linspace(-pimax, pimax, 40+1)

print('DD')
DD = two_point_counter_2d(quasar, quasar, rp_bins, pi_bins, auto_corr=True)
np.save('DD', DD.value)
print('DR')
RR = two_point_counter_2d(random, random, rp_bins, pi_bins, auto_corr=True)
np.save('RR', RR.value)
print('RR')
DR = two_point_counter_2d(quasar, random, rp_bins, pi_bins, auto_corr=False)
np.save('DR', DR.value)

