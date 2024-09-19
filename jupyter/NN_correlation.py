import numpy as np
import tqdm
from multiprocessing import Pool, Process, Pipe
from sklearn.neighbors import BallTree
from astropy.cosmology import Planck18 as cosmos
h = cosmos.H0.value/100

class two_point_counter_2d():
    def __init__(self, cata_1, cata_2, rp_bins, pi_bins, Npro=60, auto_corr=False):
        self.cata_1 = cata_1
        self.cata_2 = cata_2
        self.N1 = len(self.cata_1)
        self.N2 = len(self.cata_2)
        self.rp_bins = rp_bins
        self.pi_bins = pi_bins
        self.Npro = Npro
        self.auto_corr = auto_corr

        self.Npibins = len(pi_bins) - 1
        self.Nrpbins = len(rp_bins) - 1
        self.tree = BallTree(
            data=np.deg2rad(np.deg2rad(np.c_[cata_2['dec'], cata_2['ra']])),
            leaf_size=5, metric='haversine')          # latitude + logtitude
        self.dis_cov_1 = cosmos.comoving_distance(cata_1['z']).to('Mpc').value*h
        self.dis_cov_2 = cosmos.comoving_distance(cata_2['z']).to('Mpc').value*h
        self.weight_1 = cata_1['w']
        self.weight_2 = cata_2['w']
        self.receiver, self.sender = Pipe()
        mypool = Pool(processes=Npro)
        args = np.deg2rad(np.c_[np.deg2rad(cata_1['ra']), np.deg2rad(cata_1['dec']), self.dis_cov_1, cata_1['w']])
        assert args.shape == (len(cata_1), 4)

        progress = Process(target=self.progress_bar, args=(len(cata_1), 10))
        progress.start()
        result = mypool.map(self.stack_single_sample, args, chunksize=200)
        self.sender.send(0)
        progress.join()
        self.value = np.stack(result, axis=0)
        if auto_corr:
            self.value /= 2
    def stack_single_sample(self, arg):
        ra, dec, r, w = arg
        values = np.zeros((self.Nrpbins, self.Npibins))
        theta_bins = self.rp_bins/r
        idx_all = self.tree.query_radius(
                [[dec, ra]]*(self.Nrpbins+1), theta_bins)

        for i in range(self.Nrpbins):
            index = np.zeros(len(self.dis_cov_2), bool)
            index[idx_all[i+1]] = True
            index[idx_all[i]] = False

            dis_cov_t = self.dis_cov_2[index]                # all the points within rp[i] ~ rp[i+1]
            weight_t = self.weight_2[index]*w
            pi_bin_t = r + self.pi_bins

            pi_bin_res = np.histogram(dis_cov_t, pi_bin_t, range=(pi_bin_t[0], pi_bin_t[-1]), density=False, weights=weight_t)
            values[i,:] = pi_bin_res[0]                 # further bin the result into different pi_bins
        if self.sender is not None:
            self.sender.send(1)
        return values
    def progress_bar(self, Ntot, rate=10):
        # rate is the update rate of the progress bar.
        rate = max(1, int(rate))
        with tqdm.tqdm(total=int(Ntot/rate)) as pbar:
            count = 0
            while True:
                if self.receiver.recv():
                    count += 1
                    if count % rate == 0:
                        pbar.update(1)
                else:
                    break
            pbar.close()
def calculate_xi(DD: two_point_counter_2d, DR: two_point_counter_2d, RR: two_point_counter_2d, label):
    weight = np.fromfunction(lambda i, j: DD.weight_1[i]*DD.weight_2[j]*(i<j), shape=(DD.N1, DD.N2))
    dd = DD.value.sum(axis=0) / weight.sum()
    weight = np.fromfunction(lambda i, j: DD.weight_1[i]*DD.weight_2[j]*(i<j), shape=(RR.N1, RR.N2))
    rr = RR.value.sum(axis=0) / weight.sum()
    weight = np.fromfunction(lambda i, j: DD.weight_1[i]*DD.weight_2[j], shape=(DR.N1, DR.N2))
    dr = DR.value.sum(axis=0) / weight.sum()

    return (dd-2*dr+rr)/rr

class w_theta_counter():
    def __init__(self, autocorr, ra1, dec1, weight1, bins, ra2=None, dec2=None, weight2=None, Npro=60):
        self.autocorr = bool(autocorr)
        self.ra1 = ra1
        self.dec1 = dec1
        self.weight1 = weight1
        if not self.autocorr:
            self.ra2 = ra2
            self.dec2 = dec2
            self.weight2 = weight2
        else:
            self.ra2 = ra1
            self.dec2 = dec1
            self.weight2 = weight1
        self.tree = BallTree(\
            data=np.deg2rad(np.c_[self.dec2, self.ra2]),
            leaf_size=5, metric='haversine')          # latitude + logtitude
            
        self.bins = np.array(bins)
        self.Nbins = len(bins) - 1
        
        self.receiver, self.sender = Pipe()
        mypool = Pool(processes=Npro)
        args = np.c_[np.deg2rad(self.ra1), np.deg2rad(self.dec1), self.weight1]

        progress = Process(target=self.progress_bar, args=(len(self.ra1), 10))
        progress.start()
        result = mypool.map(self.stack_single_sample, args, chunksize=200)
        self.sender.send(0)
        progress.join()
        self.count = np.stack(result, axis=0)
        if self.autocorr:
            self.count /= 2
        if self.autocorr:
            self.normalized = self.count.sum(axis=0) / \
                (0.5*(np.sum(np.sum(self.weight1)*self.weight2)-\
                      np.sum(self.weight1[i]*self.weight2[i] for i in range(len(self.weight1)))))
        else:
            self.normalized = self.count.sum(axis=0) / \
                np.sum(np.sum(self.weight1)*self.weight2)
    def stack_single_sample(self, arg):
        ra, dec, w = arg
        values = np.zeros(self.Nbins)
        idx = self.tree.query_radius(
                [[dec, ra]]*(self.Nbins+1), self.bins)
        for i in range(self.Nbins):
            index = np.zeros(len(self.weight2), bool)
            index[idx[i+1]] = True
            index[idx[i]] = False
            weight_t = self.weight2[index].sum()*w
            values[i] = weight_t
        if self.sender is not None:
            self.sender.send(1)
        return values
    def progress_bar(self, Ntot, rate=10):
        # rate is the update rate of the progress bar.
        rate = max(1, int(rate))
        with tqdm.tqdm(total=int(Ntot/rate)) as pbar:
            count = 0
            while True:
                if self.receiver.recv():
                    count += 1
                    if count % rate == 0:
                        pbar.update(1)
                else:
                    break
            pbar.close()

if __name__ == '__main__':
    file_rand = './catalogue/random_quasar_lss.npy'
    file_data = './catalogue/quasar_lss.npy'
    data_cata = np.load(file_data)
    random_cata = np.load(file_rand)
    print('finish loading catalogue')

    sep_min = 3
    sep_max = 100
    # Nbins = 15
    rp_bins = np.geomspace(sep_min, sep_max, 15+1)
    rp = np.sqrt(rp_bins[1:]*rp_bins[:-1])
    pimax = 40
    # Npibins = 40
    pi_bins = np.linspace(-pimax, pimax, 40+1)
    t = two_point_counter_2d(data_cata, data_cata, rp_bins, pi_bins, 60, True)
    np.save('tmp', t.value)
