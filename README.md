# quasar CMB lensing correlation
Use SDSS DR16 QSO catalogue & Planck 2018 MV kappa map.

Use `healpy` to pixelize the kappa map and apply point-pixel correlation.

Check the validation of the code by comparing the result with doi:10.1093/mnras/stw2482, which uses CMASS and LOWZ, and Planck 2015 result.