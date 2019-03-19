# test likelihood notebook
%matplotlib inline
%load_ext autoreload
%autoreload 2
%config InlineBackend.figure_format = 'retina'

# %%
import numpy as np
import matplotlib.pyplot as plt
from mfpy import CrossSpectrum as CrossSpectrum

# %%
l_max = 6000
filename = 'data/wmap7_act_lcdm_bestfit_lensedCls_6000.dat'
data = np.genfromtxt(filename, unpack=True, usecols=(0,1),
                     dtype=[('ell', '<i8'),('cltt', '<f8')])
cl_th = np.zeros(data['ell'][-1]+1)
cl_th[data['ell']] = data['cltt'] # make the index match the ell
plt.plot(cl_th)
# %%



# %% initialize cross spectra
from itertools import product

root_dir = 'data/data_act/equa/'
freqs = ('148', '220')
get_freq_ind = { '148' : 0, '220' : 1 }
seasons = ('3e', '4e')

# read spectra and windows
spectra, windows, ell_lims = {}, {}, {}
cross_list = []
for f1, f2, s1, s2 in product(freqs, freqs, seasons, seasons):
    # some spectra are redundant, see Das et al 2013 Sec 4.4
    if (not (f1 == '220' and f2 == '148') and
        (not (s1 == '4e' and s2 == '3e') or f1 != f2)):
        cross_list.append( (f1, f2, s1, s2) )
        specfile = f'{root_dir}spectrum_{f1}x{f2}_season{s1}xseason{s2}.dat'
        winfile = f'{root_dir}BblMean_{f1}x{f2}_season{s1}xseason{s2}.dat'
        # store all the data from files into the dictionaries
        ell, cl, nl = np.genfromtxt(specfile, unpack=True)
        spectra[(f1, f2, s1, s2)] = CrossSpectrum(ell, cl, nl)
        windows[(f1, f2, s1, s2)] = np.genfromtxt(winfile)[:l_max+1,3:-1].T

bin_center, bin_left, bin_right = np.genfromtxt(
    f'{root_dir}binningFile.dat', unpack=True)
invcov = np.genfromtxt(f'{root_dir}Inverse_Realistic_Cov_Mat_Equa.dat')
Nbins = np.sum([len(spectra[s].ell) for s in spectra])
invcov = np.reshape(invcov, (Nbins, Nbins))
cov = np.load(f'{root_dir}Cov_Mat_Equa.npy')

# %%
import fgspectra.models

# define a list of components
components = [
    fgspectra.models.ThermalSZ(),
    fgspectra.models.KinematicSZ(),
    fgspectra.models.CIBP(),
    fgspectra.models.CIBC(),
    fgspectra.models.tSZxCIB(),
    fgspectra.models.RadioPointSources(),
    fgspectra.models.GalacticCirrus()
]

# annoying effective frequencies list: must change this interface
f0_sz     =146.9
f0_synch  =147.6
f0_dust   =149.7
f1_sz     =220.2
f1_synch  =217.6
f1_dust   =219.6
eff_freq = [[f0_sz, f1_sz], [f0_sz, f1_sz], [f0_dust, f1_dust], [f0_dust, f1_dust],
            [f0_dust, f1_dust], [f0_synch, f1_synch], [f0_dust, f1_dust] ]

# %%
test_par = {
    'a_tSZ' : 4.66,
    'a_kSZ' : 1.60,
    'a_p' : 6.87,
    'beta_p' : 2.08,
    'a_c' : 6.10,
    'beta_c' : 2.08,
    'n_CIBC' : 1.20,
    'xi' : 0.09,
    'a_s' :3.50,
    'a_g' :0.88,

    'f0_sz' :146.9,
    'f0_synch'  :147.6,
    'f0_dust'   :149.7,
    'f1_sz'     :220.2,
    'f1_synch'  :217.6,
    'f1_dust'   :219.6
}


def secondaries(par):
    """returns mixing matrix """
    mix_result = ([comp.get_mix(freqs=freq, l_max=l_max, **par)
            for comp, freq in zip(components, eff_freq)])
    mix_result = np.sum(mix_result,axis=0)
    mix_result[:,:,:2] = 0.0
    return mix_result

test_secondaries = secondaries(test_par)
# %%
for s in cross_list:
    print( windows[s].shape )
# %%
for s in cross_list:
    print( spectra[s].ell.shape )
# %%
data_vector = np.hstack( spectra[s].spectrum for s in cross_list  )
data_vector.shape

l_mins = [500] * 3 + [1500] * 7
# %%
theory_vector = np.hstack( (
    windows[s] @ ((cl_th +
    test_secondaries[
    get_freq_ind[s[0]], get_freq_ind[s[1]]
    ] ) / (np.arange(l_max+1) * (np.arange(l_max+1)+1) / 2 / np.pi + 1e-9))
)[-spectra[s].ell.shape[0]:][spectra[s].ell > l_min] for s, l_min in zip(cross_list, l_mins) )

# %%
windows[cross_list[-1]].shape
# %%
data_vector.shape
# %%

# plt.plot( theory_vector / data_vector )
# plt.plot( data_vector )
# plt.yscale("log")
# %%
plt.figure(figsize=(10,2))
plt.plot(theory_vector / data_vector)
# plt.plot( data_vector )
# %%

-2 * (data_vector - theory_vector).T @ invcov @ (data_vector - theory_vector)

# plt.plot( windows[cross_list[0]][:l_max+1,1:-4].T @ np.arange(l_max+1) )

# %%

# %%
centers.shape
