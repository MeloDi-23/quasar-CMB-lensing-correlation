#!/bin/bash
# conda activate corrfunc
python clustering.py catalogue/quasar_lss_all.npy catalogue/random_quasar_lss_all.npy 1>> auto.log
python clustering.py catalogue/quasar_lss_high_L.npy catalogue/random_quasar_lss_all.npy 1>> auto.log
python clustering.py catalogue/quasar_lss_low_L.npy catalogue/random_quasar_lss_all.npy 1>> auto.log
python clustering.py catalogue/quasar_lss_z1.npy catalogue/random_z1.npy 1>> auto.log
python clustering.py catalogue/quasar_lss_z2.npy catalogue/random_z2.npy 1>> auto.log
python clustering.py catalogue/quasar_lss_z3.npy catalogue/random_z3.npy 1>> auto.log