# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 13:47:57 2021

@author: Filip
"""

import hyperspy.api as hs
from coreshellp import CoreShellP,CoreShellSpec
import numpy as np
from specMapDiff import specMapDiff,cLoadsFacs,setCalibration,rel
from specerr import *

s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")
sAg = s.inav[0]
sCu = s.inav[-1]


size = 100
cs_mat = CoreShellP(size,30.0,20.0,1,1,1)
cs_mat.scale(1e-1)

CoreShell = CoreShellSpec(cs_mat,sCu,sAg)

core = hs.signals.Signal1D(CoreShell.core)
core = setCalibration(core,sAg)

shell = hs.signals.Signal1D(CoreShell.shell)
shell = setCalibration(shell,sAg)
#shell.add_poissonian_noise()
cs = hs.signals.Signal1D(CoreShell.getmatr())

cs = setCalibration(cs,sAg)
cs.add_poissonian_noise()
avgcounts = cs.inav[:,:].data.sum()/(cs.data.shape[0]*cs.data.shape[1])
print(avgcounts)


decomp_dim = 2
cs.decomposition(output_dimension = decomp_dim ,algorithm='NMF')
NMF_facs = cs.get_decomposition_factors()
NMF_loads = cs.get_decomposition_loadings()

cs.blind_source_separation(number_of_components=decomp_dim)#,algorithm="orthomax"
bss_facs = cs.get_bss_factors()
bss_loads = cs.get_bss_loadings()


hs.plot.plot_spectra(NMF_facs.isig[0.0:10000.0],style='cascade') 
hs.plot.plot_images(NMF_loads, cmap='mpl_colors',
            label=['Cu Core', 'Ag Shell'],
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
    
    
hs.plot.plot_spectra(bss_facs.isig[0.0:10000.0],style='cascade') 
hs.plot.plot_images(bss_loads, cmap='mpl_colors',
            axes_decor='off', per_row=1,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})

NMF = cLoadsFacs(NMF_loads,NMF_facs)
BSS = cLoadsFacs(bss_loads,bss_facs)

NMFShellDiff = specMapDiff(NMF.inav[0],shell)
NMFCoreDiff = specMapDiff(NMF.inav[1],core)
NMFShellDiff.plot()

NMFShellDiff = setCalibration(NMFShellDiff,sAg)
imNMFShell = NMFShellDiff.get_lines_intensity()

imShell = shell.get_lines_intensity()

relNMFAg = rel(imNMFShell[0],imShell[0])
relNMFCu = rel(imNMFShell[1],imShell[1])
relNMF = [relNMFAg,relNMFCu]

hs.plot.plot_images(relNMF,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    label = ['Silver Diff%','Koppar Diff%'],
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})


hs.plot.plot_images(imNMFShell,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    label = ['Silver Diff abs','Koppar Diff abs'],
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
