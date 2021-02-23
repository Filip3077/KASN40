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

# Here the two spectrums are made from simulations of pure elements and then linearly combined to the desiered 
# core and shell compositions. 
sAgPure = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")

sAg = 0.88*sAgPure.data + 0.12*sCuPure.data
sCu = 0.9*sCuPure.data + 0.1*sAgPure.data

# CoreShellP generates two 3D matrixes of a sphere. One consisting of the core and one as the shell. 
dens = 1 # The density here can be seen as having the unit nm^-3 to make the values in the sphere matrix unitless.
x = CoreShellP(50,15.0,10.0,dens,dens,1) # 50x50 pixels, 15nm outer radius, 10nm core radius, densities, 1x1 nm pixel size.

# CoreShellSpec fills the core and shell matrices from above with the simulated spectra and are then turned into 
# HyperSpy objects for the latter comparrisons. Combining these gives us a HyperSpy object of a whole particle (p). 
a = CoreShellSpec(x,sCu,sAg)
core = hs.signals.Signal1D(a.core)
shell = hs.signals.Signal1D(a.shell)
particle = core + shell

p = hs.signals.Signal1D(particle)
p.add_poissonian_noise() # Adds poissonian noise to the existing spectra.

# For nicer plots, HyperSpy needs some meta data: 
p.set_signal_type("EDS_TEM") 
p.axes_manager.signal_axes[0].units = 'keV' # Note, the unit of the energy axis is important to be able to plot.
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.add_elements(['Ag','Cu']) 
p.add_lines(['Ag_La','Cu_Ka'])

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
p.get_calibration_from(cal) # Calibration from a similar spectrum lets HyperSpy place the element lines correctly.

im = p.get_lines_intensity() # 

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%% 
# AbsCore = SpecErrAbs2D(p.data,a.core)
# AbsShell = SpecErrAbs2D(p.data,a.shell)
# print(AbsCore+AbsShell)

#%% Statistical Methods


dim = 2 # Since we only have the two elements in the particle only two components is needed for NMF and BSS.
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

#%% 
NMFspec = cLoadsFacs(loadings, factors)
NMFparticle = NMFspec.inav[0] + NMFspec.inav[1]

# errNMFtot = SpecErrAbs2D(p.data,NMFparticle.data)
errNMFcore = SpecErrAbs2D(core.data,NMFspec.inav[1].data)
errNMFshell = SpecErrAbs2D(shell.data,NMFspec.inav[0].data)
#%%

compare = specMapDiff(p,NMFparticle)


compare.set_signal_type("EDS_TEM") 
compare.axes_manager.signal_axes[0].units = 'keV' #OBS! Enheten på energiaxeln är viktig för att kunna plotta
compare.axes_manager[0].name = 'y'
compare.axes_manager[1].name = 'x'
compare.axes_manager['x'].units = 'nm'
compare.axes_manager['y'].units = 'nm'
compare.axes_manager[-1].name = 'E'
compare.add_elements(['Ag','Cu']) #Lägger in ämnen igen tydligen förs de inte med 
compare.add_lines(['Ag_La','Cu_Ka'])

compare.get_calibration_from(cal)

imComp = compare.get_lines_intensity()

comp_relAg = rel(imComp[0],im[0])
comp_relCu = rel(imComp[1],im[1])

plt.subplot(211)
comp_relAg.plot()
comp_relCu.plot()








