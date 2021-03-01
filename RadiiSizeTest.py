# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 13:29:07 2021

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
    ret = []
    result = np.zeros((len(rcPart.inav[:,0,0])))
    for i in range(len(rcPart.inav[:,0,0])):
        cName = 'Not defined'
        absErrBest = np.inf
        for j in range(len(oc)):
            absErr = SpecErrAbs2D(rcPart.inav[i],oc[j])
            if (absErr < absErrBest):
                absErrBest = absErr
                if (absErrBest < 0.5):
                    cName = oc[j].metadata.General.title
        ret.append(cName)
        result[i] =absErrBest
    return [result,ret]
    

#%% Generate particles
    
s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]
sC = hs.load("./Spectra/Carbonbackground.msa",signal_type="EDS_TEM")

size = 50
sCore = (sCu*0.9 + sAg*0.1)/100
sShell = (sCu*0.1 + sAg*+0.9)/100
carbonMat = addSpectrum(np.ones((size,size)),sC,4)
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'

cs_mat = []
oc_mat = []

amount = range(1,20)
for x in amount:
    mat = CoreShellP(size,20.0,x,1,1,1)
    prct = CoreShellSpec(mat,sCore,sShell)
    
    core = hs.signals.Signal1D(prct.core)
    core.metadata.General.title = 'Core'
    
    shell = hs.signals.Signal1D(prct.shell)
    shell.metadata.General.title = 'Shell'
    oc_mat.append([core,shell,background])
    
    
    prct = prct.getmatr() + carbonMat
    cs_mat.append(hs.signals.Signal1D(prct))
    cs_mat[-1].metadata.General.title = 'r = %d' %(x)
    
    
#%% Evalutate particles
decomp_dim = 3
save = []
result = np.empty((len(amount),decomp_dim))
ret = []
for i in range(len(cs_mat)):
    cs_mat[i].decomposition(output_dimension = decomp_dim ,algorithm='NMF')
    NMF_facs = cs_mat[i].get_decomposition_factors()
    NMF_loads = cs_mat[i].get_decomposition_loadings()
    save.append([NMF_facs,NMF_loads])
    reConst = cLoadsFacs(NMF_loads,NMF_facs)
    allt = match(reConst,oc_mat[i])
    result[i] = allt[0]
    ret.append(allt[1])

#%% PLots
x = range(1,20)
fig, axs = plt.subplots(3, sharex=True, sharey=True)
axs[0].plot(x,result[:,0])
axs[0].legend(['Loading 0'])

axs[1].plot(x,result[:,1])
axs[1].legend(['Loading 1'])

axs[2].plot(x,result[:,2])
axs[2].legend(['Loading 2'])

coreFacs = []
x = []
for i in range(len(ret)):
    for j in range(len(ret[i])):
        if ret[i][j] == 'Core':
            coreFacs.append(save[i][0].inav[j])
            x.append(i)
            break
cF = hs.stack(coreFacs)
hs.plot.plot_spectra(cF,style='cascade',padding=-1) 
plt.title('Core factors')
plt.legend([x])

