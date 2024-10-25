import numpy as np
from astropy import coordinates as coo
from astropy import units as u
import healpy as hp
import matplotlib.pyplot as plt
import multiprocessing as mp
from multiprocessing import shared_memory
from astropy.io import fits
from itertools import repeat
from functools import partial
from concurrent.futures import ThreadPoolExecutor

def weight_nan_mean(value, weight):
    masked_value = np.ma.MaskedArray(value, mask=np.isnan(value))
    return np.ma.average(masked_value, axis=0, weights=weight).data

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
"""
The function calculation jackknife variance is used for variance. It returns mean, variance and covariance.
I try to warp it in a class, but it turns out to be extremely slow to do so.
"""

def calculate_jackknife_variance_multithread(value, weight, label,
            value_r=None, weight_r=None, label_r=None,
            Npro=40, norm_cov=False):
    #global _label, _label_r, _value, _value_r, _weight, _weight_r
    
    assert value.shape == weight.shape
    assert len(label) == value.shape[0]
    if value_r is None:
        assert weight_r is None and label_r is None
        sub = False
    else:
        assert value_r.shape == weight_r.shape
        assert len(label_r) == value_r.shape[0]
        sub = True

    unique_label = np.unique(label)
    pro = mp.Pool(Npro)

    def resample(i):
        idx = label != i
        return weight_nan_mean(value[idx], weight[idx])

    def resample_sub(i):
        idx = label != i
        idx_r = label_r != i
        return weight_nan_mean(value[idx], weight[idx]) - \
            weight_nan_mean(value_r[idx_r], weight_r[idx_r])


    if sub:            # for those with random value (value = data - random)
        resample_func = lambda i_array: [resample_sub(i) for i in i_array]
        mean = weight_nan_mean(value, weight) - weight_nan_mean(value_r, weight_r)
    else:
        resample_func = lambda i_array: [resample(i) for i in i_array]
        mean = weight_nan_mean(value, weight)
    with ThreadPoolExecutor(Npro) as pool:
        futures = []
        unique_label = np.unique(label)
        bins = (np.arange(Npro+1)*len(unique_label)/Npro).astype(int)
        bins[-1] = len(unique_label)
        for i in range(Npro):
            futures.append(pool.submit(resample_func, unique_label[bins[i]:bins[i+1]]))
    res = []
    for fu in futures:
        res.extend(fu.result())

    jack_val = np.vstack(res)

    cov_mat = np.cov(jack_val, ddof=0, rowvar=False)*(len(unique_label)-1)
    std = np.sqrt([cov_mat[i,i] for i in range(cov_mat.shape[0])])
    if norm_cov:
        normalize_cov(cov_mat, std)
    return {'mean': mean, 'std': std, 'cov': cov_mat} 


def normalize_cov(cov_mat, std):
    for i in range(cov_mat.shape[0]):
        for j in range(cov_mat.shape[0]):
            cov_mat[i,j] /= std[i]*std[j]


def resample_sub_2(i):
    global _label, _value, _weight, _value_r, _label_r, _weight_r
    idx = _label != i
    idx_r = _label_r != i
    return weight_nan_mean(_value[idx], _weight[idx]) - \
        weight_nan_mean(_value_r[idx_r], _weight_r[idx_r])

def resample_2(i):
    global _label, _value, _weight
    idx = _label != i
    return weight_nan_mean(_value[idx], _weight[idx])

def calculate_jackknife_variance_global(value, weight, label,
            value_r=None, weight_r=None, label_r=None,
            Npro=40, norm_cov=False, return_jackknife=False):
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
        jack_val = np.vstack(
            pro.map(resample_sub_2, unique_label))
        mean = weight_nan_mean(value, weight) - weight_nan_mean(value_r, weight_r)
    else:
        jack_val = np.vstack(pro.map(resample_2, unique_label))
        mean = weight_nan_mean(value, weight)

    cov_mat = np.cov(jack_val, ddof=0, rowvar=False)*(len(unique_label)-1)
    std = np.sqrt([cov_mat[i,i] for i in range(cov_mat.shape[0])])
    if norm_cov:
        normalize_cov(cov_mat, std)
    return_value = {'mean': mean, 'std': std, 'cov': cov_mat}
    if return_jackknife:
        return_value['jackknife'] = jack_val
    return return_value


calculate_jackknife_variance = calculate_jackknife_variance_global