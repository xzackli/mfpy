# test likelihood notebook
%matplotlib inline
%load_ext autoreload
%autoreload 2
%config InlineBackend.figure_format = 'retina'

# %%
import numpy as np
import matplotlib.pyplot as plt
import mfl

# %%
ttmax = 6000
filename = 'data/wmap7_act_lcdm_bestfit_lensedCls_6000.dat'
data = np.genfromtxt(filename, unpack=True, usecols=(0,1),
                     dtype=[('ell', '<i8'),('cltt', '<f8')])
cl_tt = np.zeros(data['ell'][-1]+1)
cl_tt[data['ell']] = data['cltt'] # make the index match the ell

# %%
s1 = np.genfromtxt("data/data_act/equa/spectrum_148x148_season3exseason3e.dat",
    dtype=[("ell", "int"), ("cl", "float"), ('nl', 'float')])


# %%
