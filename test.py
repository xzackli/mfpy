# test likelihood notebook
%matplotlib inline
%load_ext autoreload
%autoreload 2
# %config InlineBackend.figure_format = 'retina'

# %%
import numpy as np
import matplotlib.pyplot as plt
import mfpy

# %%
ttmax = 6000
filename = 'data/wmap7_act_lcdm_bestfit_lensedCls_6000.dat'
data = np.genfromtxt(filename, unpack=True, usecols=(0,1),
                     dtype=[('ell', '<i8'),('cltt', '<f8')])
cl_tt = np.zeros(data['ell'][-1]+1)
cl_tt[data['ell']] = data['cltt'] # make the index match the ell
plt.plot(cl_tt)
# %%

class CrossSpectrum:
    def __init__(self, ell, spectrum, noise=None):
        '''Initialize the cross spectrum object.'''
        self.ell = ell
        self.spectrum = spectrum
        self.noise = noise


class MultiFrequencyLikelihood:
    def __init__(self, frequencies, seasons, spectra_dict, Bbl_dict):
        '''Initialize the object with frequencies/seasons.

        Pass in lists/tuples of frequencies and seasons. These are stored
        and used to establish the ordering in the covariance matrix, and the
        data vector. You need to read in the data yourself, and load them into
        dictionaries which correspond to the spectra.

        Parameters
        ----------
        frequencies : tuple or list of strings
            specify the frequencies used in the analysis
        seasons : tuple or list of strings
            specify the seasons used in the analysis
        spectra_dict : dictionary
            This is a dict containing tuples as keys, where the tuples are
            frequency/season combinations. i.e. for cross as an instance of
            CrossSpectrum, spectra_dict = {('148', '220', '3e', '4e') : cross}

        '''
        self.frequencies = frequencies
        self.seasons = seasons
        self.spectra_dict = spectra_dict
        self.Bbl_dict = Bbl_dict



# %% initialize cross spectra
from itertools import product

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
    windows[(f1,f2)] = np.genfromtxt(file)

# %%
likelihood = MultiFrequencyLikelihood(frequencies=freqs,
                                      seasons=seasons,
                                      spectra_dict=spectra,
                                      Bbl_dict=windows)


# %%
plt.plot( spectra[('148', '148', '4e', '4e')].ell,
spectra[('148', '148', '4e', '4e')].spectrum *
(spectra[('148', '148', '4e', '4e')].ell**4) )

# %%
windows[('148', '148')].shape
# %%


centers, lefts, rights = np.genfromtxt(f'{root_dir}binningFile.dat', unpack=True)

# %%
win = windows[('148', '148')]
# %%
plt.imshow( win[:,1:], aspect=0.001 )
plt.colorbar()
#%%
 win[:,0]

# %%
plt.plot( centers, win[:ttmax+1,1:].T @ cl_tt, label='$B_{b\ell} D_{\ell}^{\mathrm{theory}}$' )
plt.plot( spectra[('148', '148', '4e', '4e')].ell,
    (spectra[('148', '148', '4e', '4e')].spectrum *
    spectra[('148', '148', '4e', '4e')].ell**2 / 2 / np.pi), label='148x148 4e 4e' )
plt.ylabel(r'$D_{\ell}$')
plt.legend();
