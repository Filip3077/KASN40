# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:02:33 2021

@author: Lucas
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellFunctions import postStatprocess, genfullparticle
# from coreshellp import CoreShellP, CoreShellSpec
# from specerr import *
from specMapDiff import *
import numpy as np



# %% Particle generation

# Here the two spectrums are made from simulations of pure elements and then linearly combined to the desiered 
# core and shell compositions. 
sAgPure = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")

sAg = 0.9*sAgPure.data + 0.1*sCuPure.data
sCu = 0.9*sCuPure.data + 0.1*sAgPure.data


a = genfullparticle(50,20,15,sCu,sAg)

p = a.full
p.add_poissonian_noise() # Adds poissonian noise to the existing spectra.

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
p = setCalibration(p, cal)

im = p.get_lines_intensity() # 

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})



#%% Statistical Methods

dim = 3 # Since we only have the two elements in the particle only two components is needed for NMF and BSS.
p.decomposition(True,algorithm='NMF',output_dimension =dim)
factors = p.get_decomposition_factors() 
loadings =p.get_decomposition_loadings()


hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
plt.title('NMF')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_images(loadings, cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

p.blind_source_separation(number_of_components=dim)#,algorithm="orthomax"
bssfac = p.get_bss_factors()
bssload = p.get_bss_loadings()

hs.plot.plot_spectra(bssfac.isig[0.0:10000.0],style='cascade') 
plt.title('NMF+BSS with fastICA')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=8040, y=1.6, s='Cu-K$_\\alpha$', color='k')
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=1.6, s='Ag-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.text(x=930, y=1.6, s='Cu-L$_\\alpha$', color='k')

hs.plot.plot_images(bssload, cmap='mpl_colors',
                    axes_decor='off', per_row=3,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                  'left': 0.05, 'right':0.78})
#%% Postprocessing

NMFpost = postStatprocess(a.core,a.shell,factors,loadings,dim,cal)
BSSpost = postStatprocess(a.core,a.shell,bssfac,bssload,dim,cal)


