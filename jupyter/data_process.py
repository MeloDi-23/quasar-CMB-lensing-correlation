import numpy as np
from astropy import coordinates as coo
from astropy import units as u
import healpy as hp
import matplotlib.pyplot as plt
import multiprocessing as mp
from astropy.io import fits

def weight_nan_mean(value, weight):
    mean_val = np.zeros(value.shape[1])
    for i in range(len(mean_val)):
        idx = np.logical_not(np.isnan(value[:,i]))
        val = value[idx,i]
        wei = weight[idx,i]
        mean_val[i] = np.average(val, weights=wei)
    return mean_val

def load_data(file):
    if type(file) == str:
        file = open(file, 'rb')
    datas = []
    while True:
        try:
            datas.append(np.load(file))
        except EOFError:
            file.close()
            break
    return datas

def load_calculation_data(file, return_rp=False):
    """
    file is either in the form: r_p, data; or in the form: data
    """
    dat = load_data(file)

    if return_rp and len(dat) == 1:
        raise ValueError(f'There is no r_p information in file {file}')
    if return_rp:
        return dat[0:2]
    else:
        return dat[-1]

#class jackknife_sampler():
#    def __init__(self, value, weight, label,
#                value_r=None, weight_r=None, label_r=None,
#                Npro=40, norm_cov=False):
#        assert value.shape == weight.shape
#        assert len(label) == value.shape[0]
#        if value_r is None:
#            assert weight_r is None and label_r is None
#            self.sub = False
#        else:
#            assert value_r.shape == weight_r.shape
#            assert len(label_r) == value_r.shape[0]
#            self.sub = True
#
#        self.value_r = value_r
#        self.weight_r = weight_r
#        self.label_r = label_r
#        self.value = value
#        self.weight = weight
#        self.label = label
#        self.unique_label = np.unique(label)
#        pro = mp.Pool(Npro)
# 
#        if self.sub:            # for those with random value (value = data - random)
#            self.jack_val = np.vstack(pro.map(self.resample_sub, self.unique_label))
#            self.mean = weight_nan_mean(value, weight) - weight_nan_mean(value_r, weight_r)
#        else:
#            self.jack_val = np.vstack(pro.map(self.resample, self.unique_label))
#            self.mean = weight_nan_mean(value, weight)
# 
#        self.cov_mat = np.cov(self.jack_val, ddof=0, rowvar=False)*(len(self.unique_label)-1)
#        self.std = np.sqrt([self.cov_mat[i,i] for i in range(self.cov_mat.shape[0])])
#        if norm_cov:
#            self.normalize_cov()
#
#    def normalize_cov(self):
#        for i in range(self.cov_mat.shape[0]):
#            for j in range(self.cov_mat.shape[0]):
#                self.cov_mat[i,j] /= self.std[i]*self.std[j]
#    
#    def resample(self, i):
#        idx = self.label != i
#        return weight_nan_mean(self.value[idx], self.weight[idx])
#    def resample_sub(self, i):
#        idx = self.label != i
#        idx_r = self.label_r != i
#        return weight_nan_mean(self.value[idx], self.weight[idx]) - \
#            weight_nan_mean(self.value_r[idx_r], self.weight_r[idx_r])
#
"""
The function calculation jackknife variance is used for variance. It returns mean, variance and covariance.
I try to warp it in a class, but it turns out to be extremely slow to do so.
"""

def calculate_jackknife_variance(value, weight, label,
            value_r=None, weight_r=None, label_r=None,
            Npro=40, norm_cov=False):
    global _label, _label_r, _value, _value_r, _weight, _weight_r
    
    assert value.shape == weight.shape
    assert len(label) == value.shape[0]
    if value_r is None:
        assert weight_r is None and label_r is None
        sub = False
    else:
        assert value_r.shape == weight_r.shape
        assert len(label_r) == value_r.shape[0]
        sub = True

    _value_r = value_r
    _weight_r = weight_r
    _label_r = label_r
    _value = value
    _weight = weight
    _label = label
    unique_label = np.unique(label)
    pro = mp.Pool(Npro)

    if sub:            # for those with random value (value = data - random)
        jack_val = np.vstack(pro.map(resample_sub, unique_label))
        mean = weight_nan_mean(value, weight) - weight_nan_mean(value_r, weight_r)
    else:
        jack_val = np.vstack(pro.map(resample, unique_label))
        mean = weight_nan_mean(value, weight)

    cov_mat = np.cov(jack_val, ddof=0, rowvar=False)*(len(unique_label)-1)
    std = np.sqrt([cov_mat[i,i] for i in range(cov_mat.shape[0])])
    if norm_cov:
        normalize_cov(cov_mat, std)
    return {'mean': mean, 'std': std, 'cov': cov_mat} 


def normalize_cov(cov_mat, std):
    for i in range(cov_mat.shape[0]):
        for j in range(cov_mat.shape[0]):
            cov_mat[i,j] /= std[i]*std[j]


def resample(i):
    global _label, _value, _weight
    idx = _label != i
    return weight_nan_mean(_value[idx], _weight[idx])


def resample_sub(i):
    global _label, _label_r, _value, _value_r, _weight, _weight_r
    idx = _label != i
    idx_r = _label_r != i
    return weight_nan_mean(_value[idx], _weight[idx]) - \
        weight_nan_mean(_value_r[idx_r], _weight_r[idx_r])
