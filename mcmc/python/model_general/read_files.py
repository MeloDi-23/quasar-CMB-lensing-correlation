import numpy as np
from os import path
from scipy.interpolate import RegularGridInterpolator
import yaml
def read_wp(file_path):
    table_names = [
        '1h_cs', '1h_ss',               # this is N_massbin * N_rpbin
        '2h_cc', '2h_ss', '2h_cs',      # this is N_massbin * N_massbin * N_rpbin; Note that cs[i,j] = sc[j,i], cc[i,j] = cc[j,i]
        'cm', 'sm'                      # this is N_massbin * N_rpbin
    ]
    table = {k: np.load(path.join(file_path, f'{k}_wp.npy')) for k in table_names}
    table['rp'] = np.load(path.join(file_path, 'rp.npy'))
    return table

def interpolate_table(table, rp, rp0=None):
    """
    the table is calculated using rp0, and this function interpolate the result to rp.
    Because the wp is usually about power law, a linear interpolate between log rp-wp should work.
    """
    table_ret = {'rp': rp}
    if rp0 is None:
        rp0 = table['rp']
    for k, v in table.items():
        if k == 'rp':
            continue
        if v.ndim == 2:
            x = np.arange(v.shape[0])
            y = np.log(rp0)

            new_x = x
            new_y = np.log(rp)
            XX, YY = np.meshgrid(new_x, new_y, indexing='ij')
            interpolator = RegularGridInterpolator([x, y], v, bounds_error=False)
            new_v = (interpolator(np.stack((XX, YY), axis=-1)))
            assert new_v.shape == (v.shape[0], len(rp))
            table_ret[k] = new_v
            # table[k] = new_v
        elif v.ndim == 3:
            x = np.arange(v.shape[0])
            y = np.arange(v.shape[1])
            z = np.log(rp0)

            new_x = x
            new_y = y
            new_z = np.log(rp)
            XX, YY, ZZ = np.meshgrid(new_x, new_y, new_z, indexing='ij')
            interpolator = RegularGridInterpolator([x, y, z], v, bounds_error=False)
            new_v = (interpolator(np.stack((XX, YY, ZZ), axis=-1)))
            assert new_v.shape == (v.shape[0], v.shape[1], len(rp))
            table_ret[k] = new_v
            # table[k] = new_v
    return table_ret

def read_signal(file):
    ret_ls = []
    with open(file, 'rb') as f:
        while True:
            try: ret_ls.append(np.load(f))
            except EOFError:
                break
    return ret_ls

def read_cov(file):
    cov = np.load(file)
    return cov

def read_halo_mass_function(file):
    func = np.load(file)
    return func['Mass'], func['nh']

from run_mcmc import HODParameter
def read_config(file_name):
    config: list[dict[str, dict]]
    with open(file_name, 'r') as f:
        config = list(yaml.load_all(f, yaml.FullLoader))

    return HODParameter.from_config(config[0]), config[1]


import importlib.util

def read_module(module_path):
    """
    read module from path.
    """
    spec = importlib.util.spec_from_file_location("module_name", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module