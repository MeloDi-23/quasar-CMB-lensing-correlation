from emcee import EnsembleSampler
from emcee.backends import HDFBackend
import numpy as np
from read_files import *
from calc_wp import chi_2, N_c
from multiprocessing import Pool, Process
from collections import namedtuple
import numpy as np
from calc_wp import chi_2, N_c

class HODParameter:
    """
    This class handles the HOD parameter.
    It converts emcee parameter(array like) into keyword like parameter
    """
    parameter_name = [
        'lgMmin',
        'sig_lgM',
        'lgM0',
        'lgM1p',
        'alpha'
    ]
    parameter_tuple = namedtuple('parameter', parameter_name)
    def __init__(self, free_parameter=None, default_value={}, prior={}) -> None:
        self.free_parameter = np.ones(5, bool)
        if free_parameter:
            self.free_parameter[:] = free_parameter
        self.parameters = np.array(
            [ 13.500000, 0.500000, 13.500000, 14.500000, 1.000000 ]     # default value
        )
        self.parameter_priors = np.ones((5, 2), float)
        self.parameter_priors[:,0] = -np.inf
        self.parameter_priors[:,1] = np.inf             # default prior range

        if type(default_value) == dict:
            for k, v in default_value.items():
                self.set_default(k, v)
        else:
            self.parameters[:] = default_value

        for k, v in prior.items():
            self.set_prior(k, *v)
        
        for i in range(len(self.parameters)):
            if not (self.parameter_priors[i,0] <= self.parameters[i] <= self.parameter_priors[i,1]):            
                raise \
            ValueError('The default value {:.g} for parameter {} is out of prior {(:.g, :.g)}'\
                       .format(self.parameters[i], HODParameter.parameter_name[i], self.parameter_priors[i,0], self.parameter_priors[i,1]))

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
                part2 = '[{:g}-{:g}]'.format(prior[0], prior[1])
            else:
                part2 = '[fixed]'
            strings.append(part1 + part2)
        return repr.format(joiner.join(strings))
    def set_free(self, name):
        idx = HODParameter.parameter_name.index(name)
        if idx != -1:
            self.free_parameter[idx] = True

    def set_fix(self, name):
        idx = HODParameter.parameter_name.index(name)
        if idx != -1:
            self.free_parameter[idx] = False

    def set_default(self, name, value):
        idx = HODParameter.parameter_name.index(name)
        if idx != -1:
            self.parameters[idx] = value
        else:
            raise KeyError('{} not in parameters {}'.format(name, HODParameter.parameter_name))

    def set_prior(self, name, left=None, right=None):
        idx = HODParameter.parameter_name.index(name)
        if idx != -1:
            if left:
                self.parameter_priors[idx,0] = left
            if right:
                self.parameter_priors[idx,1] = right
        else:
            raise KeyError('{} not in parameters {}'.format(name, HODParameter.parameter_name))
    @property
    def ndim(self):
        return self.free_parameter.sum()
    
    def initial_parameter(self, size):
        ndim = self.ndim
        init_p = np.zeros(shape=(size, ndim), dtype=float)
        init_p0 = self.parameters[self.free_parameter]
        prior = self.parameter_priors[self.free_parameter]
        for i in range(ndim):
            sig = (prior[i,1] - prior[i,0])/6
            init_p[:,i] = init_p0[i] + np.random.uniform(-sig, sig, size)
            init_p[:,i] = np.minimum(np.maximum(init_p[:,i], prior[i,0]+1e-6), prior[i,1]-1e-6)
        return init_p

    def full_parameters(self, parameters):
        # convert from emcee para to standard para.
        ret = self.parameters.copy()
        ret[self.free_parameter] = parameters
        return HODParameter.parameter_tuple._make(ret)

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
        
        # raise NotImplementedError()
        chi2 = chi_2(signal, cov_inv, full_para, wp_table, logM, Nh)/2      # use global variables for better performance.
        if np.isnan(chi2).any():
            np.savetxt('para', full_para)
            raise Exception('chi^2 = nan.')
        return prior - chi2
"""
parameter_name = [
        'lgMmin',
        'sig_lgM',
        'lgM0',
        'lgM1p',
        'alpha'
    ]
"""

if __name__ == '__main__':
    signal_length = 12

    wp_table = read_wp('../wp_table')
    rp, signal = read_signal('../signal.npy', use_range=[0, signal_length])
    cov = read_cov('../cov.npy', use_range=[0, signal_length])
    cov_inv = np.linalg.inv(cov)
    logM, Nh = read_halo_mass_function('../halo_mass_function.npy')
    # for better performance I would recommend setting these as global variable

    sep_min = 0.1
    sep_max = 100
    sep_N = 30

    r_pbins = np.geomspace(sep_min, sep_max, sep_N+1)

    rp0 = (r_pbins[1:]*r_pbins[:-1])**0.5
    interpolate_table(wp_table, rp, rp0)
    import sys
    args = sys.argv
    if len(args) < 2:
        print('Input the config file. See config.yaml for example.')
        sys.exit(-1)
    
    parameter = HODParameter(free_parameter=[1, 1, 1, 1, 1],
                             default_value=[ 12, 1, 11.500000, 12, 1.000000 ],
                             prior={'lgMmin': (9, 16), 'sig_lgM': (0.1, 5), 'lgM0': (9, 18), 'lgM1p': (10, 16), 'alpha': (0.1, 5)})
    
    other_configs = read_config(args[1], parameter)

    print(parameter)
    Ndim = parameter.ndim
    # don't set the lgMmin too large or sig_lgM too small, otherwise N_c will be all zeros
    available_fields = ['Nwalkers', 'Nstep', 'Nburnin', 'Npro', 'backend_file', 'numpy_file']
    Nwalkers = 40
    Nstep = 4000
    Nburnin = 300
    Npro = 40
    backend_file = ''
    numpy_file = ''

    for k in available_fields:
        if k in other_configs:
            locals()[k] = other_configs[k]

    init_state = parameter.initial_parameter(Nwalkers)
    log_prob_fn = parameter.ln_prob
    pool = Pool()

    if backend_file:
        backend = HDFBackend(backend_file)
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
