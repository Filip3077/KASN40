# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 17:35:18 2021

@author: Filip
"""
import numpy as np
import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec
from specMapDiff import *

#%%
#Generation of particles
s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")
sCu = s.inav[-1]
sAg = s.inav[0]
sCore = sCu*0.9 + sAg*0.1
sShell = sAg*0.9 + sCu*0.1
totSize = 100
fullImg = np.zeros((totSize,totSize,2048))
pSize = 50
x = CoreShellP(pSize,20,16,1,1,1)
particle = CoreShellSpec(x,sCore,sShell)


for i in range(2):
    for j in range(2):
        fullImg[pSize*i:(1+i)*pSize,pSize*j:(1+j)*pSize] = particle.getmatr()
    
SI = hs.signals.Signal1D(fullImg)
SI.plot()


#%%
#NMF and other post processing  


SI = setCalibration(SI,sCu)

decomp_dim = 3
SI.decomposition(output_dimension = decomp_dim ,algorithm='NMF')
NMF_facs = SI.get_decomposition_factors()
NMF_loads = SI.get_decomposition_loadings()
 
hs.plot.plot_spectra(NMF_facs.isig[0.0:10000.0],style='cascade') 
hs.plot.plot_images(NMF_loads, cmap='mpl_colors',
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
    
SI.blind_source_separation(number_of_components=decomp_dim)#,algorithm="orthomax"
bss_facs = SI.get_bss_factors()
bss_loads = SI.get_bss_loadings()

hs.plot.plot_spectra(bss_facs,style='cascade') 
hs.plot.plot_images(bss_loads, cmap='mpl_colors',
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})