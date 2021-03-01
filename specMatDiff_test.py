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
from ErrorStack import ErrorStack
from ROI4_hyperspy import varimax
s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")
sAg = s.inav[0]
sCu = s.inav[-1]

sCore  = sCu*0.9 + sAg*0.1
sShell = sCu*0.1 + sAg*0.9

size = 100
cs_mat = CoreShellP(size,30.0,20.0,1,1,1)
cs_mat.scale(1e-1)

CoreShell = CoreShellSpec(cs_mat,sCore,sShell)

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
cs.decomposition(output_dimension = decomp_dim ,algorithm='sklearn_pca')
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

NMF = ErrorStack(NMF_facs,NMF_loads,[core,shell],['Core','Shell'])


relNMF = NMF.stack[0].relMaps + NMF.stack[1].relMaps

hs.plot.plot_images(relNMF,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})


pos = [50,50]

LC_factors = [1, 0.72980399]
Sbw = shell.inav[pos].estimate_background_windows(line_width=[5.0, 8.0])
Shellim = shell.inav[pos].get_lines_intensity(background_windows=Sbw)
shellQ = shell.inav[pos].quantification(Shellim, method='CL', factors=LC_factors,composition_units='weight')

Cbw = core.inav[pos].estimate_background_windows(line_width=[5.0, 8.0])
Ccim = core.inav[pos].get_lines_intensity(background_windows=Cbw)
coreQ = core.inav[pos].quantification(Ccim, method='CL', factors=LC_factors,composition_units='weight')


print('\n'+str(coreQ[0].data[0]) + '\n' + str(coreQ[1].data[0]) + '\n\n')
print(str(shellQ[0].data[0]) + '\n' +str(shellQ[1].data[0]))
