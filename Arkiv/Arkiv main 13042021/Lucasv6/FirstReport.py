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
'''  '''
# Here the two spectrums, made from DTSA-II simulations of pure elements, are loaded and linearly combined 
# to the desiered core and shell compositions. 
sAgPure = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")

sAg = 0.9*sAgPure.data + 0.1*sCuPure.data # Shell composition: 90%Ag-10%Cu
sCu = 0.9*sCuPure.data + 0.1*sAgPure.data # Core composition: 10%Ag-90%Cu

# CoreShellP generates two 3D matrices of a sphere. One consisting of the core and one as the shell. 
dens = 1 # The density here can be seen as having the unit nm^-3 to make the values in the sphere matrix unitless.
x = CoreShellP(50,20.0,15.0,dens,dens,1) # 50x50 pixels, 20nm outer radius, 15nm core radius, densities, 1x1 nm pixel size.

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
p.get_calibration_from(cal) # Calibration from a similar spectrum so HyperSpy place the element lines correctly.

im = p.get_lines_intensity() # 

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})



#%% Statistical Methods


dim = 2 # Since we only have the two elements in the particle only two components is needed for NMF and BSS.
p.decomposition(True,algorithm='NMF',output_dimension =dim) # The "True" variable tells the function to normalize poissonian noise.
factors = p.get_decomposition_factors() 
loadings =p.get_decomposition_loadings()

# Plotting NMF
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

p.blind_source_separation(number_of_components=dim) # BSS Based on the already performed NMF.
bssfac = p.get_bss_factors()
bssload = p.get_bss_loadings()

# Plotting BSS
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

#%% Total Sums of Errors
''' Here the total sum of errors is calculated for the different image comparisons. '''

NMFspec = cLoadsFacs(loadings, factors)
NMFparticle = NMFspec.inav[0] + NMFspec.inav[1]

errNMFtot = SpecErrAbs2D(NMFparticle,particle)
errNMFcore = SpecErrAbs2D(NMFspec.inav[0],core)
errNMFshell = SpecErrAbs2D(NMFspec.inav[1],shell)

BSSspec = cLoadsFacs(bssload,bssfac)
BSSparticle = BSSspec.inav[0] + BSSspec.inav[1]

errBSStot = SpecErrAbs2D(BSSparticle,particle)
errBSScore = SpecErrAbs2D(BSSspec.inav[0],core)
errBSSshell = SpecErrAbs2D(BSSspec.inav[1],shell)
#%% Error Maps 
''' Here error maps are calculated by subtractng the original images from the ones made from factors and loadings. '''
NMFcompareTot = specMapDiff(p,NMFparticle)
NMFcompareTot = setCalibration(NMFcompareTot, cal)
NMFcompareCore = specMapDiff(core, NMFspec.inav[0])
NMFcompareCore = setCalibration(NMFcompareCore,cal)
NMFcompareShell = specMapDiff(shell,NMFspec.inav[1])
NMFcompareShell = setCalibration(NMFcompareShell,cal)

BSScompareTot = specMapDiff(p,BSSparticle)
BSScompareTot = setCalibration(BSScompareTot, cal)
BSScompareCore = specMapDiff(core, BSSspec.inav[0])
BSScompareCore = setCalibration(BSScompareCore, cal)
BSScompareShell = specMapDiff(shell, BSSspec.inav[1])
BSScompareShell = setCalibration(BSScompareShell, cal)

#%% 
''' This just divides the maps into signals of Ag and Cu by themselfes. '''
imNMFcompTot = NMFcompareTot.get_lines_intensity()
imNMFcompCore = NMFcompareCore.get_lines_intensity()
imNMFcompShell = NMFcompareShell.get_lines_intensity()
imBSScompTot = BSScompareTot.get_lines_intensity()
imBSScompCore = BSScompareCore.get_lines_intensity()
imBSScompShell = BSScompareShell.get_lines_intensity()

#%% Relative error maps
'''  '''
particle = setCalibration(particle, cal)
imComp = particle.get_lines_intensity()

relNMFcompTot = [rel(imNMFcompTot[0], imComp[0]),rel(imNMFcompTot[1], imComp[1])]
relNMFcompCore = [rel(imNMFcompCore[0], imComp[0]),rel(imNMFcompCore[1], imComp[1])]
relNMFcompShell = [rel(imNMFcompShell[0], imComp[0]),rel(imNMFcompShell[1], imComp[1])]
relBSScompTot = [rel(imBSScompTot[0],imComp[0]),rel(imBSScompTot[1],imComp[1])]
relBSScompCore = [rel(imBSScompCore[0],imComp[0]),rel(imBSScompCore[1],imComp[1])]
relBSScompShell = [rel(imBSScompShell[0],imComp[0]),rel(imBSScompShell[1],imComp[1])]

#%% Plotting absolute maps

NMFplot = [imNMFcompTot[0], imNMFcompCore[0], imNMFcompShell[0], imNMFcompTot[1], imNMFcompCore[1], imNMFcompShell[1]]

NMFtitle = 'NMF absolute error X-ray maps'
NMFlabel = ['Ag abs.diff Tot','Ag abs.diff Core','Ag abs.diff Shell',
            'Cu abs.diff Tot','Cu abs.diff Core','Cu abs.diff Shell']
hs.plot.plot_images(NMFplot,suptitle=NMFtitle, label=NMFlabel,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

BSSplot = [imBSScompTot[0], imBSScompCore[0], imBSScompShell[0], imBSScompTot[1], imBSScompCore[1], imBSScompShell[1]]

BSStitle = 'BSS absolute error X-ray maps'
BSSlabel = ['Ag abs.diff Tot','Ag abs.diff Core','Ag abs.diff Shell',
            'Cu abs.diff Tot','Cu abs.diff Core','Cu abs.diff Shell']
hs.plot.plot_images(BSSplot,suptitle=BSStitle, label=BSSlabel,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%% Plotting relative maps

NMFplot = [relNMFcompTot[0], relNMFcompCore[0], relNMFcompShell[0], relNMFcompTot[1], relNMFcompCore[1], relNMFcompShell[1]]

NMFtitle = 'NMF relative error X-ray maps'
BSSlabel = ['Ag rel.diff Tot','Ag rel.diff Core','Ag rel.diff Shell',
            'Cu rel.diff Tot','Cu rel.diff Core','Cu rel.diff Shell']
hs.plot.plot_images(NMFplot,suptitle=NMFtitle, label=BSSlabel,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

BSSplot = [relBSScompTot[0], relBSScompCore[0], relBSScompShell[0], relBSScompTot[1], relBSScompCore[1], relBSScompShell[1]]

BSStitle = 'BSS relative error x-ray maps'
BSSlabel = ['Ag rel.diff Tot','Ag rel.diff Core','Ag rel.diff Shell',
            'Cu rel.diff Tot','Cu rel.diff Core','Cu rel.diff Shell']
hs.plot.plot_images(BSSplot,suptitle=BSStitle, label=BSSlabel,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%% Quantification
'''  '''
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

co.get_calibration_from(cal)
sh.get_calibration_from(cal)


kfactors = [0.72980399,1.0]
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





