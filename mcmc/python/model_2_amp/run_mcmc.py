from emcee import EnsembleSampler
from emcee.backends import HDFBackend
import numpy as np
from read_files import *
from calc_wp import chi_2, N_c
from multiprocessing import Pool, Process
from collections import namedtuple
from types import SimpleNamespace
import numpy as np
from calc_wp import chi_2, N_c
import os


class HODParameter:
    """
    This class handles the HOD parameter.
    It converts emcee parameter(array like) into keyword like parameter
    """
    @classmethod
    def from_config(cls, config_dict: dict[str, dict]):
        parameter_name = []
        free_parameter = []
        default_value = []
        prior = []
        for k, v in config_dict.items():
            if not isinstance(v, dict) or 'default' not in v:
                continue
            parameter_name.append(k)

            v.setdefault('free', True)
            v.setdefault('prior', [-np.inf, np.inf])

            default_value.append(v['default'])
            free_parameter.append(v['free'])
            prior.append(v['prior'])
        return cls(parameter_name, free_parameter, default_value, prior)
    def __init__(self, parameter_name=[], free_parameter=[], default_value=[], prior={}) -> None:
        """
        This code doesn't check type.
        parameter_name: list[str]
        free_parameter: list[any]
        default_value: list[float] | dict[str, float]
        prior: dict[str, list[float]] | list[list[float]]
        """
        para_dim = len(parameter_name)
        self.parameter_name = parameter_name
        self.free_parameter = np.ones(para_dim, bool)   # default: all free
        self.parameters = np.ones(para_dim)*np.nan
        self.parameter_priors = np.ones((para_dim, 2), float)
        self.parameter_priors[:,0] = -np.inf
        self.parameter_priors[:,1] = np.inf             # default prior range
        # self.parameter_tuple = namedtuple('parameter', parameter_name)
        
        if free_parameter:
            self.free_parameter[:] = free_parameter

        if type(default_value) == dict:
            for k, v in default_value.items():
                self.set_default(k, v)
        else:
            self.parameters[:] = default_value

        if type(prior) == dict:
            for k, v in prior.items():
                self.set_prior(k, *v)
        else:
            self.parameter_priors[:] = prior

        for i in range(len(self.parameters)):
            if np.isnan(self.parameters[i]):
                raise ValueError("parameter {} default value must be assigned.".format(self.parameter_name[i]))
            if not (self.parameter_priors[i,0] <= self.parameters[i] <= self.parameter_priors[i,1]):            
                raise \
            ValueError('The default value {:.g} for parameter {} is out of prior {(:.g, :.g)}'\
                       .format(self.parameters[i], self.parameter_name[i], self.parameter_priors[i,0], self.parameter_priors[i,1]))

    def __repr__(self) -> str:
        repr =      'HODParameter({})'
        joiner =  '\n             '
        strings = []
        for i in range(len(self.parameter_name)):
            name = self.parameter_name[i]
            default = self.parameters[i]
            prior = self.parameter_priors[i]
            free = self.free_parameter[i]
            part1 = '{:<10} {:<4g}'.format(name+':', default)
            if free:
                part2 = '[{:g}~{:g}]'.format(prior[0], prior[1])
            else:
                part2 = '[fixed]'
            strings.append(part1 + part2)
        return repr.format(joiner.join(strings))
    def set_free(self, name):
        idx = self.parameter_name.index(name)
        if idx != -1:
            self.free_parameter[idx] = True

    def set_fix(self, name):
        idx = self.parameter_name.index(name)
        if idx != -1:
            self.free_parameter[idx] = False

    def set_default(self, name, value):
        idx = self.parameter_name.index(name)
        if idx != -1:
            self.parameters[idx] = value
        else:
            raise KeyError('{} not in parameters {}'.format(name, self.parameter_name))

    def set_prior(self, name, left=None, right=None):
        idx = self.parameter_name.index(name)
        if idx != -1:
            if left:
                self.parameter_priors[idx,0] = left
            if right:
                self.parameter_priors[idx,1] = right
        else:
            raise KeyError('{} not in parameters {}'.format(name, self.parameter_name))
    @property
    def ndim(self):
        return self.free_parameter.sum()
    
    def initial_parameter(self, size):
        ndim = self.ndim
        init_p = np.zeros(shape=(size, ndim), dtype=float)
        init_p0 = self.parameters[self.free_parameter]
        prior = self.parameter_priors[self.free_parameter]
        for i in range(ndim):
            sig = min(prior[i,1] - init_p0[i],  init_p0[i] - prior[i,0]) / 2
            if np.isinf(sig):
                sig = 10
            init_p[:,i] = init_p0[i] + np.random.uniform(-sig, sig, size)
            # init_p[:,i] = np.minimum(np.maximum(init_p[:,i], prior[i,0]+1e-6), prior[i,1]-1e-6)
        return init_p

    def full_parameters(self, parameters):
        # convert from emcee para to standard para.
        ret = self.parameters.copy()
        ret[self.free_parameter] = parameters
        return namedtuple('par', self.parameter_name)._make(ret)

    def ln_prior(self, full_para):
        # flat prior
        if np.all((self.parameter_priors[:,0] <= full_para) & (self.parameter_priors[:,1] >= full_para)):
            if not np.isclose(N_c(logM, full_para).sum(), 0): 
                return 0
        
        return -np.inf

    def ln_prob(self, parameters):
        """
        likely hood function
        descriped as exp(-chi^2/2)
        ln prob = -chi^2/2
        """
        full_para = self.full_parameters(parameters)
        prior = self.ln_prior(full_para)
        if np.isneginf(prior):
            return prior
        
        chi2 = chi_2(signal, cov_inv, full_para, wp_table_auto, wp_table_cross, logM, Nh)/2      # use global variables for better performance.
        if np.isnan(chi2).any():
            np.savetxt('para', full_para)
            raise Exception('chi^2 = nan.')
        return prior - chi2

if __name__ == '__main__':

    import sys
    args = sys.argv
    if len(args) < 2:
        print('Input the config file. See config.yaml for example.')
        sys.exit(-1)

    with open(args[1], 'r') as f:
        par_configs, other_configs = list(yaml.load_all(f, yaml.FullLoader))
    # ========== reading configs ==========

    # don't set the lgMmin too large or sig_lgM too small, otherwise N_c will be all zeros
    available_fields = ['Nwalkers', 'Nstep', 'Nburnin', 'Npro', 'auto_range', 'cross_range',
                        'backend_file', 'numpy_file', 'wp_table_path', 'signal_path']
    Nwalkers = 40
    Nstep = 4000
    Nburnin = 300
    Npro = 40
    auto_range = None
    cross_range = None
    backend_file = ''
    numpy_file = ''
    wp_table_path = '../wp_table'
    signal_path = '../signal'

    for k in available_fields:
        if k in other_configs:
            locals()[k] = other_configs[k]

    wp_table = read_wp(wp_table_path)
    rp_auto, signal_auto, rp_cross, signal_cross, number_density = read_signal(signal_path + '/signal.npy')
    cov = read_cov(signal_path + '/cov.npy')

    logM, Nh = read_halo_mass_function('../halo_mass_function.npy')
    # for better performance I would recommend setting these as global variable

    # ========== apply the fitting range ==========

    if auto_range or cross_range:
        if auto_range is None:
            auto_range = [0, len(rp_auto)]
        if cross_range is None:
            cross_range = [0, len(rp_cross)]
        l_a, r_a = auto_range[0], auto_range[1]
        l_c, r_c = cross_range[0], cross_range[1]

        auto_size = len(rp_auto)
        cross_size = len(rp_cross)
        rp_auto = rp_auto[l_a:r_a]
        rp_cross = rp_cross[l_c:r_c]
        signal_auto = signal_auto[l_a:r_a]
        signal_cross = signal_cross[l_c:r_c]

        cov = np.vstack((
            np.hstack((cov[l_c:r_c, l_c:r_c], cov[l_c:r_c, cross_size+l_a:cross_size+r_a], cov[l_c:r_c, -1].reshape(-1, 1))),
            np.hstack((cov[cross_size+l_a:cross_size+r_a, l_c:r_c], cov[cross_size+l_a:cross_size+r_a, cross_size+l_a:cross_size+r_a], cov[cross_size+l_a:cross_size+r_a, -1].reshape(-1, 1))),
            np.hstack((cov[-1, l_c:r_c], cov[-1, cross_size+l_a:cross_size+r_a], [cov[-1, -1]]))
        ))

    signal = np.concatenate((signal_cross, signal_auto, [number_density]))
    cov_inv = np.linalg.inv(cov)

    # ========== handle interpolating table ==========
    sep_min = 0.1
    sep_max = 100
    sep_N = 30
    # sep_min = 0.9
    # sep_max = 100
    # sep_N = 20

    r_pbins = np.geomspace(sep_min, sep_max, sep_N+1)
    rp0 = (r_pbins[1:]*r_pbins[:-1])**0.5
    # the original rp for the table. (it actually equals to rp)
    wp_table_auto = interpolate_table(wp_table, rp_auto)
    wp_table_cross = interpolate_table(wp_table, rp_cross)

    parameter = HODParameter.from_config(par_configs)
    print(parameter)
    Ndim = parameter.ndim

    init_state = parameter.initial_parameter(Nwalkers)
    log_prob_fn = parameter.ln_prob
    pool = Pool()

    if backend_file:
        if os.path.exists(backend_file):
            res = input('Warning: file already exists. Overwrite it?[y/n]')
            if res not in ['y', 'Y', 'yes', 'Yes']:
                print('Exiting')
                sys.exit(1)
        backend = HDFBackend(backend_file)
        backend.reset(Nwalkers, Ndim)
        print(f'backend save to {backend_file}.')
    else:
        backend = None

    sampler = EnsembleSampler(Nwalkers, Ndim, log_prob_fn, pool=pool, backend=backend)

    sampler.run_mcmc(init_state, nsteps=Nstep, progress=True)

    if numpy_file:
        chain = sampler.get_chain(flat=True, discard=Nburnin)
        names = np.array(parameter.parameter_name)[parameter.free_parameter]
        chain_array = np.zeros(chain.shape[0], dtype=[(n, float) for n in names])
        for i, n in enumerate(names):
            chain_array[n] = chain[:,i]
        np.save(numpy_file, chain_array)
