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

#%% NMF errors
particle = core + shell
NMFspec = cLoadsFacs(loadings, factors)
NMFparticle = NMFspec.inav[0] + NMFspec.inav[1]

errNMFtot = SpecErrAbs2D(particle.data,NMFparticle.data)
errNMFcore = SpecErrAbs2D(core.data,NMFspec.inav[1].data)
errNMFshell = SpecErrAbs2D(shell.data,NMFspec.inav[0].data)
#%% BSS errors
BSSspec = cLoadsFacs(bssload,bssfac)
BSSparticle = BSSspec.inav[0] + BSSspec.inav[1]

errBSStot = SpecErrAbs2D(particle.data,BSSparticle.data)
errBSScore = SpecErrAbs2D(core.data,BSSspec.inav[0].data)
errBSSshell = SpecErrAbs2D(shell.data,BSSspec.inav[1].data)
#%%

NMFcompareTot = specMapDiff(p,NMFparticle)
NMFcompareTot = setCalibration(NMFcompareTot, cal)
NMFcompareCore = specMapDiff(core, NMFspec.inav[1])
NMFcompareCore = setCalibration(NMFcompareCore,cal)
NMFcompareShell = specMapDiff(shell,NMFspec.inav[0])
NMFcompareShell = setCalibration(NMFcompareShell,cal)

BSScompareTot = specMapDiff(p,BSSparticle)
BSScompareTot = setCalibration(BSScompareTot, cal)
BSScompareCore = specMapDiff(core, BSSspec.inav[0])
BSScompareCore = setCalibration(BSScompareCore, cal)
BSScompareShell = specMapDiff(shell, BSSspec.inav[1])
BSScompareShell = setCalibration(BSScompareShell, cal)

#%%

imNMFcompTot = NMFcompareTot.get_lines_intensity()
imNMFcompCore = NMFcompareCore.get_lines_intensity()
imNMFcompShell = NMFcompareShell.get_lines_intensity()
imBSScompTot = BSScompareTot.get_lines_intensity()
imBSScompCore = BSScompareCore.get_lines_intensity()
imBSScompShell = BSScompareShell.get_lines_intensity()

#%%
particle = setCalibration(particle, cal)
imComp = particle.get_lines_intensity()

relNMFcompTot = [rel(imNMFcompTot[0], imComp[0]),rel(imNMFcompTot[1], imComp[1])]
relNMFcompCore = [rel(imNMFcompCore[0], imComp[0]),rel(imNMFcompCore[1], imComp[1])]
relNMFcompShell = [rel(imNMFcompShell[0], imComp[0]),rel(imNMFcompShell[1], imComp[1])]
relBSScompTot = [rel(imBSScompTot[0],imComp[0]),rel(imBSScompTot[1],imComp[1])]
relBSScompCore = [rel(imBSScompCore[0],imComp[0]),rel(imBSScompCore[1],imComp[1])]
relBSScompShell = [rel(imBSScompShell[0],imComp[0]),rel(imBSScompShell[1],imComp[1])]

#%%

NMFplot = [relNMFcompTot[0], relNMFcompCore[0], relNMFcompShell[0], relNMFcompTot[1], relNMFcompCore[1], relNMFcompShell[1]]

NMFtitle = 'NMF X-ray line intensities'
BSSlabel = ['Ag rel.diff Tot','Ag rel.diff Core','Ag rel.diff Shell',
            'Cu rel.diff Tot','Cu rel.diff Core','Cu rel.diff Shell']
hs.plot.plot_images(NMFplot,suptitle=NMFtitle, label=BSSlabel,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

BSSplot = [relBSScompTot[0], relBSScompCore[0], relBSScompShell[0], relBSScompTot[1], relBSScompCore[1], relBSScompShell[1]]

BSStitle = 'BSS X-ray line intensities'
BSSlabel = ['Ag rel.diff Tot','Ag rel.diff Core','Ag rel.diff Shell',
            'Cu rel.diff Tot','Cu rel.diff Core','Cu rel.diff Shell']
hs.plot.plot_images(BSSplot,suptitle=BSStitle, label=BSSlabel,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%% Quantification

core_spec = np.sum(a.core,(0,1))
shell_spec = np.sum(a.shell,(0,1))


co = hs.signals.Signal1D(core_spec)
sh = hs.signals.Signal1D(shell_spec)
co.set_signal_type("EDS_TEM")
sh.set_signal_type("EDS_TEM")
co.axes_manager.signal_axes[0].units = 'keV'
sh.axes_manager.signal_axes[0].units = 'keV'
co.add_elements(['Ag','Cu'])
sh.add_elements(['Ag','Cu'])
co.add_lines(['Ag_La','Cu_Ka'])
sh.add_lines(['Ag_La','Cu_Ka'])

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
co.get_calibration_from(cal)
sh.get_calibration_from(cal)


kfactors = [1.0, 0.72980399]
bw = co.estimate_background_windows(line_width=[5.0, 2.0])
intensities = co.get_lines_intensity(background_windows=bw)
weight_percent_cu_core = co.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
weight_percent_ag_core = co.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

bw = sh.estimate_background_windows(line_width=[5.0, 2.0])
intensities = sh.get_lines_intensity(background_windows=bw)
weight_percent_cu_shell = sh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
weight_percent_ag_shell = sh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

NMFfac = factors.data
NMFco = hs.signals.Signal1D(NMFfac[0])
NMFsh = hs.signals.Signal1D(NMFfac[1])
NMFco.set_signal_type("EDS_TEM")
NMFsh.set_signal_type("EDS_TEM")
NMFco.axes_manager.signal_axes[0].units = 'keV'
NMFsh.axes_manager.signal_axes[0].units = 'keV'
NMFco.add_elements(['Ag','Cu'])
NMFsh.add_elements(['Ag','Cu'])
NMFco.add_lines(['Ag_La','Cu_Ka'])
NMFsh.add_lines(['Ag_La','Cu_Ka'])

NMFco.get_calibration_from(cal)
NMFsh.get_calibration_from(cal)

bw = NMFco.estimate_background_windows(line_width=[5.0, 2.0])
intensities = NMFco.get_lines_intensity(background_windows=bw)
NMFweight_percent_cu_core = NMFco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
NMFweight_percent_ag_core = NMFco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

bw = NMFsh.estimate_background_windows(line_width=[5.0, 2.0])
intensities = NMFsh.get_lines_intensity(background_windows=bw)
NMFweight_percent_cu_shell = NMFsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
NMFweight_percent_ag_shell = NMFsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

BSSfac = bssfac.data
BSSco = hs.signals.Signal1D(BSSfac[0])
BSSsh = hs.signals.Signal1D(BSSfac[1])
BSSco.set_signal_type("EDS_TEM")
BSSsh.set_signal_type("EDS_TEM")
BSSco.axes_manager.signal_axes[0].units = 'keV'
BSSsh.axes_manager.signal_axes[0].units = 'keV'
BSSco.add_elements(['Ag','Cu'])
BSSsh.add_elements(['Ag','Cu'])
BSSco.add_lines(['Ag_La','Cu_Ka'])
BSSsh.add_lines(['Ag_La','Cu_Ka'])

BSSco.get_calibration_from(cal)
BSSsh.get_calibration_from(cal)

bw = BSSco.estimate_background_windows(line_width=[5.0, 2.0])
intensities = NMFco.get_lines_intensity(background_windows=bw)
BSSweight_percent_cu_core = BSSco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
BSSweight_percent_ag_core = BSSco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

bw = BSSsh.estimate_background_windows(line_width=[5.0, 2.0])
intensities = BSSsh.get_lines_intensity(background_windows=bw)
BSSweight_percent_cu_shell = BSSsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
BSSweight_percent_ag_shell = BSSsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data





