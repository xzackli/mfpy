
# %%
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from mfpy import CrossSpectrum, MultiFrequencyLikelihood


def δ(x,y):
    return float(x == y)

def Σ(x):
    return np.sum(x)

Cbdict = {}

def Cb(α, A, β, B):
    # returns 0 if not in dict
    return Cbdict[α+A+β+B]

def Nb(α, A, β, B):
    # returns 0 if not in dict
    return Nbdict[α+A+β+B]

n_d = [4, 4, 4, 4]
nu_b = 1

def Cov_bb_Das(Cb, α, A, β, B, γ, C, τ, D):
    """Compute analytic errorbars.

    A, B, C, D correspond to frequencies,
    α, β, γ, τ correspond to seasons
    """
    N = 0.0
    term0 = 2 * Cb**2 / nu_b
    term1 = 0.0
    term2 = 0.0
    for i in range(n_d[0]):
        for j in range(n_d[1]):
            for k in range(n_d[2]):
                for l in range(n_d[3]):
                    # equation A3 in Das et al. 2013
                    N += (1 - δ(i,j) * δ(α,β)) * (1 - δ(k,l) * δ(γ,τ))

                    # equation A6 in Das et al. 2013
                    term1 += (
                        δ(i,k) * δ(A,C) * δ(α,γ) * Nb(α,A,α,A) +
                        δ(j,l) * δ(B,D) * δ(β,τ) * Nb(β,B,β,B) +
                        δ(i,l) * δ(A,D) * δ(α,τ) * Nb(α,A,α,A) +
                        δ(j,k) * δ(B,C) * δ(β,γ) * Nb(β,B,β,B)
                    ) * (1 - δ(i,j) * δ(α,β)) * (1 - δ(k,l) * δ(γ,τ))

                    term2 += Nb(α,A,α,A) * Nb(β,B,β,B) * (
                        δ(i,k) * δ(A,C) * δ(α,γ) * δ(j,l) * δ(B,D) * δ(β,τ) +
                        δ(i,l) * δ(A,D) * δ(α,τ) * δ(j,k) * δ(B,C) * δ(β,γ)
                    ) * (1 - δ(i,j) * δ(α,β)) * (1 - δ(k,l) * δ(γ,τ))

    term1 *= Cb / (N * nu_b)
    term2 *= 1 / (N * nu_b)
    print(N)
    return term0 + term1 + term2


# %%

def get_ACT_equatorial():
    root_dir = 'data/data_act/equa/'
    freqs = ('148', '220')
    seasons = ('3e', '4e')

    # read spectra
    spectra = {}
    for f1, f2, s1, s2 in product(freqs, freqs, seasons, seasons):
        if ((f1 == '220' and f2 == '148') or
            ((s1 == '4e' and s2 == '3e') and f1 == f2)):
            file = f'{root_dir}spectrum_{f2}x{f1}_season{s2}xseason{s1}.dat'
        else:
            file = f'{root_dir}spectrum_{f1}x{f2}_season{s1}xseason{s2}.dat'
        cls = np.genfromtxt(file, dtype=[('l', 'int64'), ('cl', 'd'), ('nl', 'd')])
        spectra[(f1, f2, s1, s2)] = CrossSpectrum(cls['l'], cls['cl'], cls['nl'])

    # read windows
    windows = {}
    for f1, f2 in product(('148', '220'), repeat=2):
        file = f'{root_dir}BblMean_{f1}x{f2}_season4exseason4e.dat'
        if f1 == '220' and f2 == '148':
            file = f'{root_dir}BblMean_{f2}x{f1}_season4exseason4e.dat'
        windows[(f1,f2)] = np.genfromtxt(file)[:,1:].T


    likelihood = MultiFrequencyLikelihood(frequencies=freqs,
                                          seasons=seasons,
                                          spectra_dict=spectra,
                                          Bbl_dict=windows)
    return likelihood
# %%

highell_likelihood_dict = {
    'amp_tsz' :4.66e0,  #tSZ power
    'amp_ksz' :1.60e0,  #kSZ power
    'xi'      :0.09e0,  #tSZ-CIB cross-correlation coefficient
    'amp_d'   :6.87e0,  #CIB Poisson power
    'amp_s'   :3.50e0,  #Radio Poisson power in ACT
    'amp_s2'  :1.37e0,  #Radio Poisson power in SPT
    'amp_c'   :6.10e0,  #CIB clustered power
    'n'       :1.20e0,  #CIB spectral index in ell
    'beta_c'  :2.08e0,  #CIB clustered emissivity index (for T:9.7K)
    'beta_d'  :2.08e0,  #CIB Poisson emissivity index (for T:9.7K)
    'alpha_s' :-0.5e0,  #Radio spectral index in flux
    'amp_gs'  :0.31e0,  #ACT-S Cirrus power
    'amp_ge'  :0.88e0,  #ACT-E Cirrus power

    'cas1'    :1.01e0,  #ACT-S 148 GHz calibration
    'cas2'    :1.02e0,  #ACT-S 218 GHz calibration
    'cae1'    :1.00e0,  #ACT-E 148 GHz calibration
    'cae2'    :0.99e0,  #ACT-E 218 GHz calibration
    'cal_1'   :1.00e0,  #SPT 95 GHz calibration
    'cal_2'   :1.01e0,  #SPT 150 GHz calibration
    'cal_3'   :1.02e0   #SPT 220 GHz calibration
}
