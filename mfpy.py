from itertools import product
import numpy as np
import fgspectra.models

class CrossSpectrum:
    def __init__(self, ell, spectrum, l_min=2, l_max=10000, noise=None):
        '''Initialize the cross spectrum object.'''
        self.ell = ell
        self.spectrum = spectrum
        self.noise = noise
        self.l_min = l_min
        self.l_max = l_max


class WindowFunction:
    def __init__(self, ell, window):
        '''Initialize the cross spectrum object.'''
        self.ell = ell
        self.window = window

def pad_arr_ell2(arr):
    """
    Adds two columns of zeros to a 2D array on the last dimension,
    to fix arrays which start at ell=2 so that they start at ell=0 (the
    convention used in this code.)
    """
    result_size = np.array(arr.shape)
    result_size[-1] += 2
    result = np.zeros(result_size)
    result[...,2:] = arr
    return result

class ACT_E:
    def __init__(self, th_l_max=6000, root_dir = 'data/data_act/equa/'):
        '''Initialize the object with the ACT-E Files.

        Parameters
        ----------
        l_max : int
            maximum multipole to compute likelihoods for
        root_dir : string
            where to look for files
        '''
        # %% initialize cross spectra
        self.freqs = ('148', '220')
        self.get_freq_ind = { '148' : 0, '220' : 1 }
        self.seasons = ('3e', '4e')
        self.l_max = 10000

        # read spectra and windows
        self.spectra = {}
        self.windows = {}
        self.cross_list = []
        Nbins = 0 # counter for number of spectra bins

        ell_keep_stack = []

        for f1, f2, s1, s2 in product(self.freqs, self.freqs, self.seasons, self.seasons):
            # some spectra are redundant, see Das et al 2013 Sec 4.4
            if (not (f1 == '220' and f2 == '148') and
                (not (s1 == '4e' and s2 == '3e') or f1 != f2)):
                self.cross_list.append( (f1, f2, s1, s2) )
                specfile = f'{root_dir}spectrum_{f1}x{f2}_season{s1}xseason{s2}.dat'
                winfile = f'{root_dir}BblMean_{f1}x{f2}_season4exseason4e.dat'
                # winfile = f'{root_dir}BblMean_{f1}x{f2}_season{s1}xseason{s2}.dat'

                # store all the data from files into the dictionaries
                spec_l_min = 500 if (f1=='148' and f2=='148') else 1500
                spec_l_max = 10000
                ell, cl, nl = np.genfromtxt(specfile, unpack=True)
                ell_keep_stack.append(np.logical_and(ell>spec_l_min, ell<th_l_max))
                self.spectra[(f1, f2, s1, s2)] = CrossSpectrum(
                    ell, cl,
                    l_min=spec_l_min, l_max=spec_l_max, noise=nl)
                Nbins += len(ell)

                # first column of win_data is ells, rest is matrix
                win_data = np.genfromtxt(winfile).T
                self.windows[(f1, f2, s1, s2)] = WindowFunction(
                    ell=np.arange(self.l_max+1),
                    window=pad_arr_ell2(win_data[1:,:self.l_max-1]))

        self.good_bins = np.hstack(ell_keep_stack)
        self.bin_left, self.bin_right, self.bin_center = np.genfromtxt(
            f'{root_dir}binningFile.dat', unpack=True)
        invcov = np.genfromtxt(f'{root_dir}Inverse_Realistic_Cov_Mat_Equa.dat')
        self.invcov = np.reshape(invcov, (Nbins, Nbins))
        self.cov = np.load(f'{root_dir}Cov_Mat_Equa.npy')

        # define a list of components
        self.fg_components = [
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
        self.eff_freq = [[f0_sz, f1_sz], [f0_sz, f1_sz], [f0_dust, f1_dust], [f0_dust, f1_dust],
                    [f0_dust, f1_dust], [f0_synch, f1_synch], [f0_dust, f1_dust] ]

        self.eff_par = {
            'effective_frequencies' : True,
            'f0_sz' : f0_sz, 'f0_synch'  :f0_synch, 'f0_dust'   :f0_dust,
            'f1_sz'     :f1_sz, 'f1_synch'  :f1_synch, 'f1_dust'   :f1_dust
        }

        # CALIBRATIONS
        self.cal = [1.0, 0.99]

    def get_sec(self, par):
        """Return a mixing matrix with the secondary components.

        This is for frequencies 148x148, 148x220, and 220x220.
        """
        par_ = par.copy()
        par_.update(self.eff_par)  # put effective freqs is in there

        mix_result = ([comp.get_mix(freqs=freq, l_max=self.l_max, **par_)
            for comp, freq in zip(self.fg_components, self.eff_freq)])
        mix_result = np.sum(mix_result,axis=0)
        mix_result[:,:,:2] = 0.0 # zero out the ell=0, ell=1
        return mix_result


    def get_data_vector(self, l_max):
        """Stack the data Cls into one big vector for likelihood estimation."""
        data_vec = []
        for s in self.cross_list:
            spec = self.spectra[s]
            bincut = np.logical_and(spec.ell > spec.l_min,
                spec.ell < min(l_max,spec.l_max))
            data_vec.append(spec.spectrum[bincut])
        return np.hstack(data_vec)

    def get_data_ell_vector(self):
        """Stack the ells of all the spectra into one big vector."""
        return np.hstack(self.spectra[s].ell for s in self.cross_list)

    def get_theory_vector(self, cltt, sec_par, l_max, bin=True):
        """Assemble a theory Cl vector, given input spectrum and FG params."""
        sec = self.get_sec(sec_par)
        cl_list = []

        l_max = min(l_max, self.l_max)
        dl_factor = np.arange(l_max+1) * np.arange(l_max+1) / 2 / np.pi
        dl_factor[0] = 1
        for s in self.cross_list:
            if bin:
                bincut = np.logical_and(self.bin_center > self.spectra[s].l_min,
                    self.bin_center < l_max)
                i = self.get_freq_ind[s[0]]
                j = self.get_freq_ind[s[1]]
                theory_spec_binned = (
                    (self.windows[s].window[:,:l_max+1] @
                    ((cltt + sec[i,j,:l_max+1]) / dl_factor))[bincut])
                theory_spec_binned /= self.cal[i] * self.cal[j]
                cl_list.append(theory_spec_binned)
            else:
                i = self.get_freq_ind[s[0]]
                j = self.get_freq_ind[s[1]]
                cl_list.append((cltt + sec[i,j,:l_max+1])/ dl_factor)

        return np.hstack(cl_list)
