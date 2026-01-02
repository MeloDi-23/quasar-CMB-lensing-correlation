# !/bin/bash
# Run all stacking analyses and log output to run.log
python kappa_stacking_BallTree.py catalogue/quasar_lss_high_L.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log 
python kappa_stacking_BallTree.py catalogue/quasar_lss_low_L.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log
python kappa_stacking_BallTree.py catalogue/quasar_lss_all.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log
python kappa_stacking_BallTree.py catalogue/quasar_lss_z1.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log 
python kappa_stacking_BallTree.py catalogue/quasar_lss_z2.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log
python kappa_stacking_BallTree.py catalogue/quasar_lss_z3.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log
python kappa_stacking_BallTree.py catalogue/random_quasar_lss_all.npy catalogue/CMB_smoothed_6.npy -m3 -M100 -N15 1>> run.log
# shuffle test
python kappa_stacking_BallTree.py catalogue/quasar_lss_all.npy catalogue/CMB_smoothed_6_shuffle.npy -m3 -M100 -N15 1>> run.log
python kappa_stacking_BallTree.py catalogue/random_quasar_lss_all.npy catalogue/CMB_smoothed_6_shuffle.npy -m3 -M100 -N15 1>> run.log
