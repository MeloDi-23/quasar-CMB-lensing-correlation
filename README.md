# quasar CMB lensing correlation
SDSS DR16 QSO catalogue & Planck 2018 MV kappa map. Use `healpy` to pixelize the kappa map and apply point-pixel correlation.

## settings
cosmology: follows BigMDPL setting, uses
```python
from astropy.cosmology import FlatLambdaCDM
cosmos = FlatLambdaCDM(H0=67.77, Om0=0.307115, Ob0=0.048206, Tcmb0=2.7255)
```
current setting: sep = 3 ~ 100, N = 15


## current catalogue

catalogue/quasar and random catalogue, mainly used is 
- quasar_lss_[high,low]_L.npy
- quasar_lss_all.npy
- quasar_lss_z[1-3].npy
- random_quasar_lss_all.npy
- random_z[1-3].npy

produced by jupyter/catalogue_process/read_quasar_lss.ipynb and jupyter/catalogue_process/split_quasar_[L,z].ipynb

catalogue/kappa map, mainly used is
- CMB_smoothed_6.npy                      smooth kernel 6 arcmin.
- CMB_smoothed_6_shuffle.npy              shuffled kappa for testing.

produced by jupyter/catalogue_process/read_planck.ipynb

## current pipeline:

- use clustering.py for auto correlation calculation (save to auto_correlation/)
- use kappa_stacking_BallTree.py for quasar-kappa cross correlation calculation (save to calculation_data/)
- use jupyter/prepare_for_mcmc_small_scale.ipynb to perform scale cut, add small scale, recalculate the cov matrix, 
output to mcmc/python/signal_[...] dirs
- use mcmc/python/model_general/run_mcmc.py to do MCMC. the HOD setting (e.g. shape of Ns Nc function) is in HOD dir; the program takes in config files in configs dir. 

delta/gaus/step means the shape of Nc; pow/fsat means the shape of Ns/Nc (power law or constant). It is not very proper to use z1~z3 for mcmc, because the redshift range is different.

## MCMC table
MCMC interpolation table is generated using make_table_rearange, a renewed version for make_table. make_cluster_*.py does the basic number counts, and wp*.ipynb output the final table. bin/ dir is middle output from the py script, and wp/ dir is the final output dir.

the mainly used config is config_small_step_pow.yaml, use the full result. the cut on auto correlation is needed, because the last few points of auto-corr is not proper to use.

## figures
The figures used in the paper are all in figure/ dir, and the program files in figureplot/ dir.
