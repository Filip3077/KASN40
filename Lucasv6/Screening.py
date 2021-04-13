# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 11:12:05 2021

@author: Lucas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from coreshellfunctions import *

#%%

sCu = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM") #s.inav[-1]
sAg  = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM") #s.inav[0]
sC = hs.load("../Spectra/Carbonbackground.msa",signal_type="EDS_TEM")
cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")

#%%

# Design Matrix
D = [[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]

# I = Intensity (counts/pixel); C = Composition; R = Radii
# LI = 1/45; LC = 60/40; LR = 5/20
# HI = 1/2; HC = 90/10; HR = 15/20

# sShell = 0.9*sAg.data + 0.1*sCu.data
# sCore = 0.9*sCu.data + 0.1*sAg.data

# Intensity = 1/2
# a = genfullparticle(50,20,15,sCore,sShell,True, Intensity, Intensity)
# a.add_background(sC, Intensity)
# p = a.core + a.shell
# p.add_poissonian_noise()
# p = setCalibration(p, cal)
# p = p.isig[2000.0::]

#%%

size = 50
im = []; cores = []; shells = []
for i in range(8):
    #Low
    if D[i][0] == -1:
        I = 1/45
    if D[i][1] == -1:
        C = [0.6,0.4]
    if D[i][2] == -1:
        R = 5
    #High
    if D[i][0] == 1:
        I = 1/2
    if D[i][1] == 1:
        C = [0.9,0.1]
    if D[i][2] == 1:
        R = 15
    specCore = sCu*C[0]+sAg*C[1]
    specShell = sAg*C[0]+sCu*C[1]
    a = genfullparticle(size, 20, R, specCore, specShell,True, I)
    cores.append(setCalibration(a.core, cal).isig[2000.0::]); shells.append(setCalibration(a.shell, cal).isig[2000.0::])
    a.add_background(sC, I)
    p = a.core + a.shell
    p.add_poissonian_noise()
    p = setCalibration(p, cal)
    im.append(p.isig[2000.0::])

#%%

nfac = 2
facs = []; loads = []
for i in range(8):
    im[i].decomposition(True, algorithm = 'NMF', output_dimension=nfac)
    facs.append(im[i].get_decomposition_factors())
    loads.append(im[i].get_decomposition_loadings())
    
#%%

facstack = hs.stack(facs)
hs.plot.plot_spectra(facstack,style='cascade')
plt.title('NMF factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
# plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
# plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_images(loads,suptitle='NMF Loadings', cmap='mpl_colors',
                    axes_decor='off', per_row=nfac,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                              'left': 0.05, 'right':0.78})

#%%

from loadassign import checkLoadFit

kfac = [1,0.72980399]
wtCores = []; wtShells = []
for i in range(8):
    c, s = checkLoadFit(cores[i], shells[i], facs[i], loads[i])
    wtCores.append(quantify(facs[i].inav[c],kfac,[],cal.isig[2000.0::]))
    wtShells.append(quantify(facs[i].inav[s],kfac,[],cal.isig[2000.0::]))



