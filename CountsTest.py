# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 15:11:31 2021

@author: Filip
"""

#%% Some functions

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

def match(rcPart, oc):
    splitrcPart = rcPart.split()
    ret = []
    result = np.zeros((len(splitrcPart)))
    for i in range(len(splitrcPart)):
        cName = 'Not defined'
        absErrBest = np.inf
        for j in range(len(oc)):
            absErr = SpecErrAbs2D(splitrcPart[i],oc[j])
            if (absErr < absErrBest):
                absErrBest = absErr
                cName = oc[j].metadata.General.title
        ret.append(cName)
        result[i] =absErrBest
    return [result,ret]
    
    

#%% Generate base particle
s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]
sC = hs.load("./Spectra/Carbonbackground.msa",signal_type="EDS_TEM")

size = 50
sCore = (sCu*0.9 + sAg*0.1)/100
sShell = (sCu*0.1 + sAg*+0.9)/100
carbonMat = addSpectrum(np.ones((size,size)),sC,1)
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'

mat = CoreShellP(size,20.0,15.0,1,1,1)
prct = CoreShellSpec(mat,sCore,sShell)


core = hs.signals.Signal1D(prct.core)
core.metadata.General.title = 'Core'
shell = hs.signals.Signal1D(prct.shell)
shell.metadata.General.title = 'Shell'
oc = [core,shell,background]

prct = prct.getmatr() #+ carbonMat
base = hs.signals.Signal1D(prct)
base.add_poissonian_noise(keep_dtype=True)

avgcounts = base.inav[:,:].data.sum()/(base.data.shape[0]*base.data.shape[1])
base = base/(avgcounts/10)
avgcounts = base.inav[:,:].data.sum()/(base.data.shape[0]*base.data.shape[1])
print("Medelantal counts: "+str(avgcounts))


base.set_signal_type("EDS_TEM")
base.get_calibration_from(sAg)
base.add_elements(['Ag','Cu','C']) #Lägger in element igen tydligen förs de inte med 
base.add_lines(['Ag_La','Cu_Ka','C_Ka'])
im = base.get_lines_intensity()
hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})    

#%% Generate all particles 
    
amount = range(1,1100,100)     
#base = base.isig[500::]
cs = []    
decomp_dim = 3
save = []
result = []
ret = []
for fac in amount:
    cs.append(base*fac)
    cs[-1].decomposition(output_dimension = decomp_dim ,algorithm='NMF')
    NMF_facs = cs[-1].get_decomposition_factors()
    NMF_loads = cs[-1].get_decomposition_loadings()
    save.append([NMF_facs,NMF_loads])
    reConst = cLoadsFacs(NMF_loads,NMF_facs)
    #result.append(match(reConst,oc))
    

l = []
for d in save:
    temp = d[1].split()
    l += temp
    
hs.plot.plot_images(l, cmap='mpl_colors',
            axes_decor='off', per_row=3,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
    

k = []
for d in save:
    temp = d[0].split()
    k += temp
    
hs.plot.plot_spectra(k,style='cascade',padding=-1)