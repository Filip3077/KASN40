# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 14:42:45 2021

@author: Lucas
"""

import coreshellFunctions
from specMapDiff import *
import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np

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

def match(Statspec,oc):
    core = oc[0]
    shell = oc[1]
    background = oc[2]
    statcoretest = []
    statshelltest = []
    statcore = 1000
    statshell = 1000
    
    for i in range(components):
        statcoretest.append(SpecErrAbs2D(Statspec.inav[i], core))
        statshelltest.append(SpecErrAbs2D(Statspec.inav[i], shell))
        if statcoretest[i] < statcore:
            statcore = statcoretest[i]
            c = i
        if statshelltest[i] < statshell:
            statshell = statshelltest[i]
            s = i

#%% Particle generation

s = hs.load("../Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")
sCu = s.inav[-1]
sAg  = s.inav[0]
sC = hs.load("../Spectra/Carbonbackground.msa",signal_type="EDS_TEM")

size = 50
sCore = (sCu*0.9 + sAg*0.1)/10
sShell = (sCu*0.1 + sAg*+0.9)/10
carbonMat = addSpectrum(np.ones((size,size)),sC,4)
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'

cs_mat = []
oc_mat = []

for x in range(1,20,2):
    prct = genfullparticle(size,20,x,sCore,sShell)
    core = prct.core
    core.metadata.General.title = 'Core'
    shell = prct.shell
    shell.metadata.General.title = 'Shell'
    oc_mat.append([core,shell,background])
    
    prct = prct.full + carbonMat
    cs_mat.append(hs.signals.Signal1D(prct))
    cs_mat[-1].metadata.General.title = 'r = %d' %(x)

#%% Evalutate particles

decomp_dim = 3
save = []
result = np.empty((len(amount),decomp_dim))
ret = []
for i in range(len(cs_mat)):
    cs_mat[i].decomposition(True)#output_dimension = decomp_dim ,algorithm='NMF')
    facs = cs_mat[i].get_decomposition_factors()
    facs_select = facs.inav[0:decomp_dim]
    loads = cs_mat[i].get_decomposition_loadings()
    loads_select = loads.inav[0:decomp_dim]
    [rotated_facs, rotated_loads] = FullVarimax(facs_select, loads_select)#,decomp_dim)
    save.append([rotated_facs,rotated_loads])
    reConst = cLoadsFacs(rotated_loads,rotated_facs)
    allt = match(reConst,oc_mat[i])
    result[i] = allt[0]
    ret.append(allt[1])
    
#%%
    