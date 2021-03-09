# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:51:41 2021

@author: Filip
"""

import hyperspy.api as hs
import numpy as np
from coreshellp import *
from specerr import SpecErrAbs2D
from specMapDiff import cLoadsFacs
import matplotlib.pyplot as plt

def addSpectrum(a,spec,specthickness):
    L = len(spec.data)
    mat = np.zeros((len(a),len(a[0]),L))

    for i in range(0,len(a)):
        for j in range(0,len(a[0])):
            mat[i,j,0:L]=(a[i,j]/specthickness)*spec.data
            
    specMat = hs.signals.Signal1D(mat)
    specMat.set_signal_type("EDS_TEM")
    specMat.get_calibration_from(spec)
    return specMat

s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]
sC = hs.load("./Spectra/Carbonbackground.msa",signal_type="EDS_TEM")

size = 100
sCore = (sCu*0.9 + sAg*0.1)/10
sShell = (sCu*0.1 + sAg*+0.9)/10

carbonMat = addSpectrum(np.ones((size,size)),sC,1)
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'

mat = CoreShellP(size,30,20,1,1,1)
cs_mat = CoreShellSpec(mat,sCore,sShell)
cs_mat = cs_mat.getmatr() + carbonMat
amount = 10

hsList = []


for x in range(1,amount):
    sList = []
    for y in range(0,x):
        sList.append(hs.signals.Signal1D(cs_mat))
        sList[-1].add_poissonian_noise(keep_dtype=True)
    if x > 1:
        hsList.append(hs.stack(sList))
    else:
        hsList.append(sList[-1])

NMFresults = []
decomp_dim = 3
for stack in hsList:
    stack.decomposition(output_dimension = decomp_dim ,algorithm='NMF',normalize_poissonian_noise=True)
    NMF_facs = stack.get_decomposition_factors()
    NMF_loads = stack.get_decomposition_loadings()
    NMFresults.append([NMF_facs,NMF_loads])
    
#%%

for each in NMFresults:
    NMF_facs = each[0]    
    for f in NMF_facs:
        f.data /= f.data.max()
    hs.plot.plot_spectra(NMF_facs,style='cascade')
    hs.plot.plot_images(each[1].split(), cmap='mpl_colors',
            axes_decor='off', per_row=3,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
