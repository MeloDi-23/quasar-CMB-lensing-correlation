import numpy as np
import tqdm
from multiprocessing import Pool, Process, Pipe
from sklearn.neighbors import BallTree
file_rand = './catalogue/random_sample_quasar_narrow_z_cut.npy'
file_data = './catalogue/quasar_narrow_z_cut.npy'

data_cata = np.load(file_data)
random_cata = np.load(file_rand)
sep_min = 0.00092039
sep_max = 0.03067955
Nbins = 15
theta_bins = np.geomspace(sep_min, sep_max, Nbins+1)
theta = np.sqrt(theta_bins[1:]*theta_bins[:-1])

tree_r = BallTree(
    data=np.deg2rad(np.c_[random_cata['dec'], random_cata['ra']]),
    leaf_size=5, metric='haversine')          # latitude + logtitude
tree_d = BallTree(
    data=np.deg2rad(np.c_[data_cata['dec'], data_cata['ra']]),
    leaf_size=5, metric='haversine')          # latitude + logtitude
Npro = 60
cata = random_cata
tree = tree_r
output = './calculation_data/output_rr'


def stack_single_sample(arg):
    global sender, tree
    r"""
    To stack around a single position
    Estimator
    \Sigma_l(r) = (\sum_p w_lp \kappa_p \Sigma_c)/(\sum_p w_lp)
    w_lp = w_l * Sigma_c^(-2)
    \Sigma_c = c^c/(4\pi G) \chi_s/(\chi_l (\chi_s-\chi_l) (1+z_l))
    """
    ra, dec = arg[0], arg[1]
    values = np.zeros(Nbins)

    count_all = tree.query_radius(
            [[dec, ra]]*(Nbins+1), theta_bins, count_only=True)
    for i in range(Nbins):
        count_in, count_out = count_all[i], count_all[i+1]
        count = count_out - count_in
        values[i] = count
    if sender is not None:
        sender.send(1)

    return values


def progress_bar(receiver, Ntot, rate=10):
    # rate is the update rate of the progress bar.
    rate = max(1, int(rate))
    with tqdm.tqdm(total=int(Ntot/rate)) as pbar:
        count = 0
        while True:
            if receiver.recv():
                count += 1
                if count % rate == 0:
                    pbar.update(1)
            else:
                break
        pbar.close()


receiver, sender = Pipe()
mypool = Pool(processes=Npro)


args = np.deg2rad(np.c_[cata['ra'], cata['dec']])
assert args.shape == (len(cata), 2)

progress = Process(target=progress_bar, args=(receiver, len(cata), 10))
progress.start()
result = mypool.map(stack_single_sample, args, chunksize=200)
sender.send(0)
progress.join()
value = np.vstack(result)
np.save(output, value)
