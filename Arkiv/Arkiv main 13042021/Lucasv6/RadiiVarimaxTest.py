# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 14:20:26 2021

@author: Lucas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from coreshellfunctions import *
from specerr import *

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

def match(reConst,oc):
    ret = []
    result = np.zeros((len(reConst.inav[:,0,0])))
    for i in range(len(result)):
        cName = 'Not Defined'
        absErrBest = np.inf
        for j in range(len(oc)):
            # reConst.inav[i] = reConst.inav[i]/np.max(reConst.inav[i].data)
            # oc[j] = oc[j]/np.max(oc[j].data)
            absErr = SpecErrAbs2D(reConst.inav[i],oc[j])
            # print(oc[j].metadata.General.title+': '+str(absErr))
            if (absErr < absErrBest):
                absErrBest = absErr
                # print('New best: '+oc[j].metadata.General.title+': '+str(absErrBest))
                if absErrBest < 1.0:
                    cName = oc[j].metadata.General.title
        # print(str(i+1)+' best: '+cName+' '+str(absErrBest))
        ret.append(cName)
        result[i] = absErrBest
    return result, ret
    
# Quantification
def quantify(factors, kfac, calSpec):
    ''' Quantifies the elements in the core or shell using choosen k-factors. '''
    factors = setCalibration(factors, calSpec)

    bw = factors.estimate_background_windows(line_width=[5.0, 7.0])
    intensities = factors.get_lines_intensity(background_windows=bw)
    wtCu = factors.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data
    wtAg = factors.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data
    return wtCu, wtAg

#%% Particle generation
# s = hs.load("../Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM") #s.inav[-1]
sAg  = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM") #s.inav[0]
sC = hs.load("../Spectra/Carbonbackground.msa",signal_type="EDS_TEM")
cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")

size = 50
sCore = (sCu.data*0.9 + sAg.data*0.1)
sShell = (sCu.data*0.1 + sAg.data*+0.9)
carbonMat = addSpectrum(np.ones((size,size)),sC,4)
# cs_background = hs.signals.Signal1D(carbonMat)
# oc_background = cs_background.metadata.General.title = 'Background'
background = hs.signals.Signal1D(carbonMat)
background.metadata.General.title = 'Background'

cs_mat = []
oc_mat = []

amount = range(2,18,1)
for x in amount:
    p = genfullparticle(size, 20, x, sCore, sShell)
    core = p.core
    core.metadata.General.title = 'Core'
    shell = p.shell
    shell.metadata.General.title = 'Shell'
    oc_mat.append([core,shell,background])
    
    prct = p.full + background
    prct.add_poissonian_noise()
    # cs_mat.append(hs.signals.Signal1D(prct))
    cs_mat.append(prct)
    cs_mat[-1] = setCalibration(cs_mat[-1], cal) # Denna är visst väldigt viktig för att få Varimax att funka
    cs_mat[-1].metadata.General.title = 'r = %d' %(x)
    
#%% Varimax

decomp_dim = 2
facload = []
result = np.empty((len(amount),decomp_dim))
ret = []

for i in range(0,len(cs_mat)):
    cs_mat[i].decomposition(True)
    facs = cs_mat[i].get_decomposition_factors().inav[0:decomp_dim]
    loads = cs_mat[i].get_decomposition_loadings().inav[0:decomp_dim]
    facs_rot, loads_rot = FullVarimax(facs, loads, decomp_dim)
    facload.append([facs_rot,loads_rot])
    reConst = cLoadsFacs(loads_rot, facs_rot)
    result[i], cName = match(reConst,oc_mat[i])
    ret.append(cName)


    
#%%

# keV = np.arange(-0.2,20.28,0.01)

# for i in range(len(cs_mat)):
#     fig1, ax1 = plt.subplots()
#     for j in range(decomp_dim):
#         ax1.plot(keV,facload[i][0].inav[j].data,label=('C'+str(i+1)))
#         ax1.legend()
#     plt.title('Rotated factors for: '+cs_mat[i].metadata.General.title)
    
#     fig2, tup = plt.subplots(nrows=1,ncols=decomp_dim, figsize=(20,4))
#     for j in range(decomp_dim):
#         tup[j].imshow(facload[i][1].inav[j].data, cmap=plt.get_cmap('viridis'),vmin=-2,vmax=7)
    
# hs.plot.plot_spectra(facs_rot.isig[0.0:10000.0],style='cascade')
# plt.title('Varimax Factors')
# plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
# plt.axvline(8040, c='k', ls=':', lw=0.5)
# plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
# plt.axvline(930, c='k', ls=':', lw=0.5)
# plt.axvline(2984, c='k', ls=':', lw=0.5)
# plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')
   
# hs.plot.plot_images(loads_rot,suptitle='Varimax Loadings', cmap='mpl_colors',
#                     axes_decor='off', per_row=3,
#                     scalebar=[0], scalebar_color='white',
#                     padding={'top': 0.95, 'bottom': 0.05,
#                               'left': 0.05, 'right':0.78})

#%% Lucas Plotts

fig1, axs = plt.subplots(nrows=decomp_dim, ncols=1, sharex=True, sharey=True)
x = amount
for i in range(decomp_dim):
    axs[i].plot(x,result[:,i])#,label=('Loading '+str(i+1)))
    axs[i].title.set_text('Error in factor '+str(i+1))
plt.xlabel('Core Radii')

kfacs = [1,0.72980399]

CoreQ = []
ShellQ = []
coreF = []
shellF = []
backF = []

for i in range(len(ret)):
    
    for j in range(len(ret[i])):
        if ret[i][j] == 'Core':
            coreF.append(facload[i][0].inav[j])
            # CoreQ.append(quantify(facload[i][0].inav[j], kfacs, cal))
        elif ret[i][j] == 'Shell':
            shellF.append(facload[i][0].inav[j])
            # ShellQ.append(quantify(facload[i][0].inav[j], kfacs, cal))
        elif ret[i][j] == 'Background':
            backF.append(facload[i][0].inav[j])

            
cF = hs.stack(coreF)
sF = hs.stack(shellF)
# bF = hs.stack(backF)

hs.plot.plot_spectra(cF.isig[0.0:10000.0],style='cascade')
plt.title('Core Factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_spectra(sF.isig[0.0:10000.0],style='cascade')
plt.title('Shell Factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

# hs.plot.plot_spectra(bF.isig[0.0:10000.0],style='cascade')
# plt.title('Background Factors')
# plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
# plt.axvline(8040, c='k', ls=':', lw=0.5)
# plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
# plt.axvline(930, c='k', ls=':', lw=0.5)
# plt.axvline(2984, c='k', ls=':', lw=0.5)
# plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

        

#%% Filip plotts

x = amount
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
            coreFacs.append(facload[i][0].inav[j])
            absErr.append(round(result[i,j],2)*100)
            lText.append(str(i)+' AbsErr = '+str(absErr))
            x2.append(i)
        elif ret[i][j] == 'Shell':
            shellFacs.append(facload[i][0].inav[j])
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
