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

def checkLoadFit(core,shell,statfac,statload,components=2):
    
   '''Takes the core and shell of a core-shell particle (Hyperspy-signals) as well as \n
   the factors and loadings of a statistical method and returns the indices\n
   of the factor+loading combo that best correspond to core and shell respectively.
    '''
   Statspec = cLoadsFacs(statload, statfac)
   statcoretest = []
   statshelltest = []
   statcore = 1000
   statshell = 1000
        
   for i in range(components):
       statcoretest=SpecErrAbs2D(Statspec.inav[i], core)
       statshelltest=SpecErrAbs2D(Statspec.inav[i], shell)
       if statcoretest < statcore:
        statcore = statcoretest
        index_c = i
       if statshelltest < statshell:
            statshell = statshelltest
            index_s = i
   if index_c == index_s:
           print('Warning! It guessed the same component for core and shell.')
   return index_c,index_s #index för den som bäst passar core respektive shell

#%% Generate particles
s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]
sC = hs.load("./Spectra/Carbonbackground.msa",signal_type="EDS_TEM")

size = 100
sCore = (sCu*0.9 + sAg*0.1)/12
sShell = (sCu*0.1 + sAg*+0.9)/12
carbonMat = addSpectrum(np.ones((size,size)),sC,1)
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'


outer =  np.linspace(5,50,10)
InnerIncrement = 5

cs_mat = []
oc_mat = []


for s in outer:
    inner = np.linspace(s/InnerIncrement,s,InnerIncrement)
    innerList = []
    innerOC = []
    for c in inner:
        mat = CoreShellP(size,s,c,1,1,1)
        prct = CoreShellSpec(mat,sCore,sShell)
        core = hs.signals.Signal1D(prct.core)
        core.metadata.General.title = 'Core [S:%d C:%d]' %(s,c)
        
        shell = hs.signals.Signal1D(prct.shell)
        shell.metadata.General.title = 'Shell[S:%d C:%d]' %(s,c)
        innerOC.append([core,shell])
        
        prct = prct.getmatr() + carbonMat
        innerList.append(hs.signals.Signal1D(prct))
        innerList[-1].metadata.General.title = 'radii = %d , %d' %(s,c)
        innerList[-1].add_poissonian_noise(keep_dtype=True)
        
    oc_mat.append(innerOC)
    cs_mat.append(innerList)
        
#amount = range(1,20,2)
#for x in amount:
#    mat = CoreShellP(size,20.0,x,1,1,1)
#    prct = CoreShellSpec(mat,sCore,sShell)
#    
#    core = hs.signals.Signal1D(prct.core)
#    core.metadata.General.title = 'Core'
#    
#    shell = hs.signals.Signal1D(prct.shell)
#    shell.metadata.General.title = 'Shell'
#    oc_mat.append([core,shell,background])
#    
#    
#    prct = prct.getmatr() + carbonMat
#    cs_mat.append(hs.signals.Signal1D(prct))
#    cs_mat[-1].metadata.General.title = 'r = %d' %(x)
#    cs_mat[-1].add_poissonian_noise(keep_dtype=True)
#    
    
#%% Evalutate particles
decomp_dim = 3
save = []
result = np.empty((len(cs_mat),len(cs_mat[0]),decomp_dim))
ret = []

for s in range(len(cs_mat)+1):
    for c in range(len(cs_mat[s])):
        cs_mat[s][c].decomposition(output_dimension = decomp_dim ,algorithm='NMF',normalize_poissonian_noise=True)
        cs_mat[s][c].blind_source_separation(number_of_components=decomp_dim)
        NMF_facs = cs_mat[s][c].get_bss_factors()
        NMF_loads = cs_mat[s][c].get_bss_loadings()
        save.append([NMF_facs,NMF_loads])
        reConst = cLoadsFacs(NMF_loads,NMF_facs)
        allt = match(reConst,oc_mat[s][c])
        result[s][c] = allt[0]
        ret.append(allt[1])

#%% PLots
    
x = np.linspace(5,50,10)
y = np.linspace(1/5,1,5)

K, I =np.meshgrid(y,x)



plt.figure()
ax = plt.contour(K,I,result[:,:,0])
ax.title = 'Loading 0'
cbar = plt.colorbar(ax)

plt.figure()
ax1 = plt.contour(K,I,result[:,:,1])
ax1.title = 'Loading 1'
cbar = plt.colorbar(ax1)      

plt.figure()
ax2 = plt.contour(K,I,result[:,:,2])
ax2.title = 'Loading 2'
cbar = plt.colorbar(ax2)      

l = []
for s in save:
    l.append(s[1])

for x in range(5,50,5):
    r = l[x-5:x]
    hs.plot.plot_images(r, cmap='mpl_colors',
            axes_decor='off', per_row=3,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
        
        
        
        
        
        
        
#test = cs_mat[5]
#test.set_signal_type("EDS_TEM")
#test.get_calibration_from(sAg)
#test.add_elements(['Ag','Cu','C']) #Lägger in element igen tydligen förs de inte med 
#test.add_lines(['Ag_La','Cu_Ka','C_Ka'])
#im = test.get_lines_intensity()
#hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
#    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
#    scalebar_color='black', suptitle_fontsize=16,
#    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
#             'right':0.85, 'wspace':0.20, 'hspace':0.10})
#
#
#avgcounts = test.inav[:,:].data.sum()/(test.data.shape[0]*test.data.shape[1])
#print("Medelantal counts: "+str(avgcounts))
#    
#l = []
#for d in save:
#    temp = d[1].split()
#    l += temp
#    
#hs.plot.plot_images(l, cmap='mpl_colors',
#            axes_decor='off', per_row=3,
#            scalebar=[0], scalebar_color='white',
#            padding={'top': 0.95, 'bottom': 0.05,
#                     'left': 0.05, 'right':0.78})
#
#x = range(1,20,2)
#fig, axs = plt.subplots(3, sharex=True, sharey=True)
#axs[0].plot(x,result[:,0])
#axs[0].legend(['Loading 0'])
#
#axs[1].plot(x,result[:,1])
#axs[1].legend(['Loading 1'])
#
#axs[2].plot(x,result[:,2])
#axs[2].legend(['Loading 2'])
#
#coreFacs = []
#lText = []
#x2 = []
#absErr =[]
#
#shellFacs = []
#sText = []
#x3 = []
#absErrshell = []
#
#for i in range(len(ret)):
#    for j in range(len(ret[i])):
#        if ret[i][j] == 'Core':
#            coreFacs.append(save[i][0].inav[j])
#            absErr.append(round(result[i,j],2)*100)
#            lText.append(str(i)+' AbsErr = '+str(absErr[-1]))
#            x2.append(i)
#        elif ret[i][j] == 'Shell':
#            shellFacs.append(save[i][0].inav[j])
#            absErrshell.append(round(result[i,j],2)*100)
#            sText.append(str(i)+' AbsErr = '+str(absErrshell[-1]))
#            x3.append(i)
#
#kfacs = [1,0.72980399]
#            
#cF = hs.stack(coreFacs)
#cF.set_signal_type("EDS_TEM")
#cF.get_calibration_from(sAg)
#cF.add_elements(['Cu','Ag'])
#cF.add_lines()
#bg = cF.estimate_background_windows(line_width=[5.0, 7.0])
#intensities = cF.get_lines_intensity(background_windows=bg)
#CoreQ = cF.quantification(intensities, method='CL', factors=kfacs,composition_units='weight')
#
#
#sF = hs.stack(shellFacs)
#sF.set_signal_type("EDS_TEM")
#sF.get_calibration_from(sAg)
#sF.add_elements(['Cu','Ag'])
#sF.add_lines()
#
#sbg = sF.estimate_background_windows(line_width=[5.0, 7.0])
#intensities = sF.get_lines_intensity(background_windows=sbg)
#ShellQ = sF.quantification(intensities, method='CL', factors=kfacs,composition_units='weight')
#
#
#fig = plt.figure()
#ax = plt.subplot(111)
#hs.plot.plot_spectra(cF,style='cascade',padding=-1,fig=fig,ax=ax) 
#fig.title ='Core factors'
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(lText,loc='center left', bbox_to_anchor=(1, 0.5))
#
#
#fig = plt.figure()
#ax = plt.subplot(111)
#hs.plot.plot_spectra(sF,style='cascade',padding=-1,fig=fig,ax=ax) 
#fig.title ='Shell factors'
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(sText,loc='center left', bbox_to_anchor=(1, 0.5))
#
#
#x2 = np.array(x2)
#x3 = np.array(x3)
#
#fig2 = plt.figure()
#ax2 = plt.subplot(211)
#ax2.plot(x2,CoreQ[0].data,x2,CoreQ[1].data,x2,absErr,'k')
#ax2.legend(['wt% Ag','wt% Cu','Abs Error (%)'])
#
#ax3 = plt.subplot(212)
#ax3.plot(x3,ShellQ[0].data,x3,ShellQ[1].data,x3,absErrshell,'k')
#ax3.legend(['wt% Ag','wt% Cu','Abs Error (%)'])
