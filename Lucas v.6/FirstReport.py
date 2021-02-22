# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:21:16 2021

@author: Filip Hallböök & Lucas Nobrant 
"""
import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec
from specerr import *
from specMapDiff import *
import numpy as np



# %% Particle generation

sAgPure = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
sCuPure = hs.load("../Spectra/20nm cube Cu80Ag20.msa",signal_type="EDS_TEM")

sAg = sAgPure #0.88*sAgPure.data + 0.12*sCuPure.data
sCu = sCuPure #0.9*sCuPure.data + 0.1*sAgPure.data

dens = 1
x = CoreShellP(50,20.0,15.0,dens,dens,1)

a = CoreShellSpec(x,sCu,sAg)
core = hs.signals.Signal1D(a.core)
shell = hs.signals.Signal1D(a.shell)
particle = core + shell

p = hs.signals.Signal1D(particle)
p.add_poissonian_noise()
p.set_signal_type("EDS_TEM") 
p.axes_manager.signal_axes[0].units = 'keV' #OBS! Enheten på energiaxeln är viktig för att kunna plotta
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.add_elements(['Ag','Cu']) #Lägger in ämnen igen tydligen förs de inte med 
p.add_lines(['Ag_La','Cu_Ka'])

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
p.get_calibration_from(cal)

im = p.get_lines_intensity()

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%% 
AbsCore = SpecErrAbs2D(p.data,a.core)
AbsShell = SpecErrAbs2D(p.data,a.shell)
print(AbsCore+AbsShell)

#%% 
dim = 2
p.decomposition(True,algorithm='NMF',output_dimension =dim)
factors = p.get_decomposition_factors() #loadings är det återskapade bilderna baserat på faktorerna 7
# for f in factors:

#     f.data /= f.data.max()
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

#%% 
NMFspec = cLoadsFacs(loadings, factors)
NMFparticle = NMFspec.inav[0] + NMFspec.inav[1]

absNMF = SpecErrAbs2D(p.data,NMFparticle.data)
#%%

comp = specMapDiff(p,NMFparticle)


comp.set_signal_type("EDS_TEM") 
comp.axes_manager.signal_axes[0].units = 'keV' #OBS! Enheten på energiaxeln är viktig för att kunna plotta
comp.axes_manager[0].name = 'y'
comp.axes_manager[1].name = 'x'
comp.axes_manager['x'].units = 'nm'
comp.axes_manager['y'].units = 'nm'
comp.axes_manager[-1].name = 'E'
comp.add_elements(['Ag','Cu']) #Lägger in ämnen igen tydligen förs de inte med 
comp.add_lines(['Ag_La','Cu_Ka'])

comp.get_calibration_from(cal)

imComp = comp.get_lines_intensity()

comp_relAg = rel(imComp[0],im[0])
comp_relCu = rel(imComp[1],im[1])

plt.subplot(211)
comp_relAg.plot()
comp_relCu.plot()








