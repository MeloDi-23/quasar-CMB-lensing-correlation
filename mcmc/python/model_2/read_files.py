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
    return {k: np.load(path.join(file_path, f'{k}_wp.npy')) for k in table_names}

def interpolate_table(table, rp, rp0):
    """
    the table is calculated using rp0, and this function interpolate the result to rp.
    Because the wp is usually about power law, a linear interpolate between log rp-wp should work.
    """
    for k, v in table.items():
        if v.ndim == 2:
            x = np.arange(v.shape[0])
            y = np.log(rp0)

            new_x = x
            new_y = np.log(rp)
            XX, YY = np.meshgrid(new_x, new_y, indexing='ij')
            interpolator = RegularGridInterpolator([x, y], v, bounds_error=False)
            new_v = (interpolator(np.stack((XX, YY), axis=-1)))
            assert new_v.shape == (v.shape[0], len(rp))
            table[k] = new_v
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
            table[k] = new_v

def read_signal(file):
    signal = np.load(file)
    return signal[0], np.concatenate([signal[1], signal[2]])

def read_cov(file):
    cov = np.load(file)
    return cov

def read_halo_mass_function(file):
    func = np.load(file)
    return func['Mass'], func['nh']

from run_mcmc import HODParameter
def read_config(file_name, parameter: HODParameter):
    config: list[dict[str, dict]]
    with open(file_name, 'r') as f:
        config = list(yaml.load_all(f, yaml.FullLoader))
    for k, v in config[0].items():
        try:
            v.setdefault('free', False)
            if v['free']:
                parameter.set_free(k)
            else:
                parameter.set_fix(k)
            if 'default' in v:
                parameter.set_default(k, v['default'])
            if 'prior' in v:
                parameter.set_prior(k, *v['prior'])
        except KeyError as e:
            print(e)
    return config[1]