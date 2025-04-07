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
        
        chi2 = chi_2(signal, cov_inv, full_para, wp_table, logM, Nh)/2      # use global variables for better performance.
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
    
    # parameter = HODParameter(
    #                     parameter_name=['lgMmin', 'sig_lgM', 'lgM0', 'lgM1p', 'alpha', 'Amp'],
    #                     free_parameter=[1, 1, 1, 1, 1, 1],
    #                          default_value=[ 12, 1, 11.500000, 12, 1.000000, 0.7 ],
    #                          prior={'lgMmin': (9, 16), 'sig_lgM': (0.1, 5), 'lgM0': (9, 18), 'lgM1p': (10, 16), 'alpha': (0.1, 5), 'Amp': (0, 1)})
    # print(parameter)
    with open(args[1], 'r') as f:
        par_configs, other_configs = list(yaml.load_all(f, yaml.FullLoader))
    # ========== reading configs ==========

    # don't set the lgMmin too large or sig_lgM too small, otherwise N_c will be all zeros
    available_fields = ['Nwalkers', 'Nstep', 'Nburnin', 'Npro', 'fitting_range',
                        'backend_file', 'numpy_file', 'wp_table_path', 'signal_path']
    Nwalkers = 40
    Nstep = 4000
    Nburnin = 300
    Npro = 40
    fitting_range = None
    backend_file = ''
    numpy_file = ''
    wp_table_path = '../wp_table'
    signal_path = '../signal'

    for k in available_fields:
        if k in other_configs:
            locals()[k] = other_configs[k]

    wp_table = read_wp(wp_table_path)
    rp, signal = read_signal(signal_path + '/signal.npy')
    cov = read_cov(signal_path + '/cov.npy')
    logM, Nh = read_halo_mass_function('../halo_mass_function.npy')
    # for better performance I would recommend setting these as global variable

    # ========== apply the fitting range ==========
    if fitting_range:
        left, right = fitting_range[0], fitting_range[1]
        signal_size = len(rp)
        rp = rp[left:right]
        signal = np.hstack((signal[left:right], signal[signal_size+left:signal_size+right]))
        cov = np.vstack((
            np.hstack((cov[left:right, left:right], cov[left:right, signal_size+left:signal_size+right])),
            np.hstack((cov[signal_size+left:signal_size+right, left:right], cov[signal_size+left:signal_size+right, signal_size+left:signal_size+right]))
        ))
    cov_inv = np.linalg.inv(cov)

    # ========== handle interpolating table ==========
    sep_min = 0.1
    sep_max = 100
    sep_N = 30
    r_pbins = np.geomspace(sep_min, sep_max, sep_N+1)
    rp0 = (r_pbins[1:]*r_pbins[:-1])**0.5
    # the original rp for the table. (it actually equals to rp)
    interpolate_table(wp_table, rp, rp0)

    parameter = HODParameter.from_config(par_configs)
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
