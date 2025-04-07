from scipy.special import erf
import numpy as np
def N_c(logM, par):
    return par.Amp*0.5*(1 + erf((logM - par.lgMmin)/par.sig_lgM))

def N_s_to_N_c(logM, par):
    # N_s(M) / N_c(M)
    res = np.zeros_like(logM)
    valid = logM >= par.lgM0
    res[valid] = ((10**logM[valid] - 10**par.lgM0)/(10**par.lgM1p))**par.alpha
    return res

def N_s(logM, par):
    return N_s_to_N_c(logM, par)*N_c(logM, par)

def N_g(logM, nh, par):
    Nc = N_c(logM, par)
    Ns = Nc * N_s_to_N_c(logM, par)
    return np.sum(nh*(Nc + Ns))

def N_g_sat(logM, nh, par):
    return np.sum(nh*N_s(logM, par))


def w_p_1h(nh, Nc, Ns, ng, wp_cs, wp_ss):
    cs = np.sum(np.einsum('i,ij', nh/(ng*ng)*Nc*Ns, wp_cs))
    ss = np.sum(np.einsum('i,ij', nh/(ng*ng)*Ns*(Ns-1), wp_ss))
    return cs*2 + ss

def w_p_2h(nh, Nc, Ns, ng, wp_cc, wp_cs, wp_ss):
    cc = np.einsum('i, j, ijk', nh*Nc, nh*Nc, wp_cc)/ng**2
    cs = np.einsum('i,j,ijk', nh*Nc, nh*Ns, wp_cs)/ng**2
    ss = np.einsum('i,j,ijk', nh*Ns, nh*Ns, wp_ss)/ng**2
    return cc + cs*2 + ss

def w_p_matter(nh, Nc, Ns, ng, wp_cm, wp_sm):
    cm = np.einsum('i,ij', nh/ng*Nc, wp_cm)
    sm = np.einsum('i,ij', nh/ng*Ns, wp_sm)
    return cm+sm

def w_p(logM, nh, par, wp_table):
    Nc = N_c(logM, par)
    Ns = Nc * N_s_to_N_c(logM, par)
    ng = np.sum(nh*(Nc + Ns))
    # expression for auto correlation (g-g)
    w_p_std = w_p_1h(nh, Nc, Ns, ng, wp_table['1h_cs'], wp_table['1h_ss']) + \
        w_p_2h(nh, Nc, Ns, ng, wp_table['2h_cc'], wp_table['2h_cs'], wp_table['2h_ss'])
    # expression for cross correlation (g-m)
    w_p_m = w_p_matter(nh, Nc, Ns, ng, wp_table['cm'], wp_table['sm'])
    return np.concatenate((w_p_m, w_p_std))

def chi_2(signal, cov_inv, par, wp_table, logM, nh):
    delta = signal - w_p(logM, nh, par, wp_table)

    return np.einsum('i,ij,j', delta, cov_inv, delta)
