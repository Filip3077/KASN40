# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:50:02 2021

@author: Filip
"""

import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec
from addSpectrum import addSpectrum
import numpy as np

s = hs.load("../Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-2]
sAg  = s.inav[1]
sCarbon = hs.load("../Spectra/Carbonbackground.msa",signal_type="EDS_TEM")

size = 200
cs_mat = CoreShellP(size,20.0e-9,15.0e-9,10.49,8.96,1e-9)

csShell = addSpectrum(cs_mat.shell,sAg,1e-6)
csShell.add_poissonian_noise(keep_dtype=True)

csCore = addSpectrum(cs_mat.core,sCu,1e-6)
csCore.add_poissonian_noise(keep_dtype=True)

cMat = np.zeros((size,size,len(sCarbon.data)))

sCarbon.add_poissonian_noise(keep_dtype=True)

cs = csCore+csShell + sCarbon

cs.plot()

cs.axes_manager[0].name = 'y'
cs.axes_manager[1].name = 'x'
cs.axes_manager['x'].units = 'nm'
cs.axes_manager['y'].units = 'nm'
cs.axes_manager[-1].name = 'E'
cs.add_elements(['Ag','Cu','C']) #Lägger in element igen tydligen förs de inte med 
cs.add_lines(['Ag_La','Cu_Ka','C_Ka'])

im = cs.get_lines_intensity()
hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})