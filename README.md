# quasar CMB lensing correlation
Use SDSS DR16 QSO catalogue & Planck 2018 MV kappa map.

Use `healpy` to pixelize the kappa map and apply point-pixel correlation.

Check the validation of the code by comparing the result with doi:10.1093/mnras/stw2482, which uses CMASS and LOWZ, and Planck 2015 result.

The main code is in `kappa_stacking_BallTree.py`, in which "BallTree" means the algorithm used to find point-pixel pairs. Other code is specialized for certain purposes.

There is some jupyter notebook that generates the catalogue and the pixelized kappa map we need.

jupyter starting with "result" is used for processing the calculation result, which uses the code from `data_process.py`. I am currently trying to appply different cut to the quasar sample and see if there is any correlation on it.


## The test done:
has random
- quasar all 2048
- quasar lss 2048 1024, and l<400, 1024

no random:
- quasar z>2 2048
- quasar lowz (z<1) (2048+nbins=15, 8)
- quasar lss all l < 4096, but use a kernel to smooth

compare:
- quasar all vs quasar cut                        quasar sample doesn't matter
- quasar lss 2048 vs 1024, 2048 vs 400            cmb kappa map doesn't matter
- quasar high z low z
- quasar luminosity cut
- quasar volumne limited
