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
from specMapDiff import *
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
                if (absErrBest < 1.0):
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
sCore = (sCu*0.9 + sAg*0.1)/10
sShell = (sCu*0.1 + sAg*+0.9)/10
carbonMat = addSpectrum(np.ones((size,size)),sC,4)
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'

cs_mat = []
oc_mat = []

amount = range(1,20,2)
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

#%% PLots
x = range(1,20,2)
fig, axs = plt.subplots(3, sharex=True, sharey=True)
axs[0].plot(x,result[:,0])
axs[0].legend(['Loading 0'])

axs[1].plot(x,result[:,1])
axs[1].legend(['Loading 1'])

axs[2].plot(x,result[:,2])
axs[2].legend(['Loading 2'])

coreFacs = []
lText = []
x2 = []
absErr =[]

shellFacs = []
sText = []
x3 = []
absErrshell = []

for i in range(len(ret)):
    for j in range(len(ret[i])):
        if ret[i][j] == 'Core':
            coreFacs.append(save[i][0].inav[j])
            absErr.append(round(result[i,j],2)*100)
            lText.append(str(i)+' AbsErr = '+str(absErr))
            x2.append(i)
        elif ret[i][j] == 'Shell':
            shellFacs.append(save[i][0].inav[j])
            absErrshell.append(round(result[i,j],2)*100)
            lText.append(str(i)+' AbsErr = '+str(absErrshell))
            x3.append(i)

kfacs = [1,0.72980399]
            
cF = hs.stack(coreFacs)
cF.set_signal_type("EDS_TEM")
cF.get_calibration_from(sAg)
cF.add_elements(['Cu','Ag'])
cF.add_lines()
bg = cF.estimate_background_windows(line_width=[5.0, 7.0])
intensities = cF.get_lines_intensity(background_windows=bg)
CoreQ = cF.quantification(intensities, method='CL', factors=kfacs,composition_units='weight')


sF = hs.stack(shellFacs)
sF.set_signal_type("EDS_TEM")
sF.get_calibration_from(sAg)
sF.add_elements(['Cu','Ag'])
sF.add_lines()

sbg = sF.estimate_background_windows(line_width=[5.0, 7.0])
intensities = sF.get_lines_intensity(background_windows=sbg)
ShellQ = sF.quantification(intensities, method='CL', factors=kfacs,composition_units='weight')


fig = plt.figure()
ax = plt.subplot(111)
hs.plot.plot_spectra(cF,style='cascade',padding=-1,fig=fig,ax=ax) 
fig.title ='Core factors'
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(lText,loc='center left', bbox_to_anchor=(1, 0.5))


x2 = np.array(x2)
x3 = np.array(x3)

fig2 = plt.figure()
ax2 = plt.subplot(211)
ax2.plot(x2,CoreQ[0].data,x2,CoreQ[1].data,x2,absErr,'k')
ax2.legend(['wt% Ag','wt% Cu','Abs Error (%)'])

ax3 = plt.subplot(212)
ax3.plot(x3,ShellQ[0].data,x3,ShellQ[1].data,x3,absErrshell,'k')
ax3.legend(['wt% Ag','wt% Cu','Abs Error (%)'])
