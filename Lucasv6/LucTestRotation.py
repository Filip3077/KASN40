# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:05:58 2021

@author: Lucas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from coreshellfunctions import *

#%% Particle generation

core_spect = hs.load('../Spectra/core.msa').data
shell_spect = hs.load('../Spectra/shell.msa').data

size = 50
a = genfullparticle(size, 20, 15, core_spect, shell_spect)
s = a.full
s.add_poissonian_noise()

cal = hs.load("../Spectra/shell.msa",signal_type="EDS_TEM")
s = setCalibration(s, cal)

s.change_dtype('float64')

ld = np.abs(np.linspace(0,size-1,size).reshape((size,1))@np.ones(size).reshape((1,size))-(size-1)/2)

d_core = np.ones((50,50))*15
d_shell = 5

l_core = np.sqrt(d_core - np.sqrt(ld**2 + ld.T**2))
l_core = np.nan_to_num(l_core)

l_shell = np.sqrt((d_core+d_shell )- np.sqrt(ld**2 + ld.T**2))-l_core
l_shell = np.nan_to_num(l_shell)

c = np.ones((size,2,2048))
c[:,0,:] = core_spect
c[:,1,:] = shell_spect

l=np.transpose(np.array([l_core, l_shell]))*0.5

# plt.imshow(l[:,:,1])

data_true = l@c  # Vet inte riktigt vad l@c inebär...
data = np.random.poisson(data_true)
s2 = hs.signals.Signal1D(data)

s2 = setCalibration(s2, cal)
s2.change_dtype('float64')

#%% PCA

s.decomposition(True)
s.plot_explained_variance_ratio(n=20, vline=True)

s2.decomposition(True)
s2.plot_explained_variance_ratio(n=20, vline=True)

nfac = 2 #antalet faktorer att ta med

#Välj nfac antal komponenter att gå vidare med
s1_factors = s.get_decomposition_factors().inav[0:nfac]
s1_loadings = s.get_decomposition_loadings().inav[0:nfac]
s2_factors = s2.get_decomposition_factors().inav[0:nfac]
s2_loadings = s2.get_decomposition_loadings().inav[0:nfac]
allFacs = hs.signals.Signal1D(np.random.random((4, 2048)),signal_type="EDS_TEM")
allFacs.inav[0:nfac] = s1_factors; allFacs.inav[nfac:nfac+2] = s2_factors
allFacs = setCalibration(allFacs, cal)

# Kör våra vanliga plottar
hs.plot.plot_spectra(allFacs,style='cascade')
plt.title('SVD factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')


hs.plot.plot_images([s1_loadings, s2_loadings],suptitle='SVD Loadings', cmap='mpl_colors',
                    axes_decor='off', per_row=2,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                              'left': 0.05, 'right':0.78})

#%% Minimizing "compactness"

'''
Jag förutsätter här att själva metoden för att hitta kompakthet stämmer även om jag inte hittat något online som 
direkt motsvarar detta.
'''
from scipy.optimize import minimize_scalar
s1_loads = (s1_loadings.data[0,:,:], s1_loadings.data[1,:,:])
s1_facs = (s1_factors.data[0,:], s1_factors.data[1,:])
s2_loads = (s2_loadings.data[0,:,:], s2_loadings.data[1,:,:])
s2_facs = (s2_factors.data[0,:], s2_factors.data[1,:])

def RotCompactness(a, *loads):

    return np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1]).sum(axis=(0,1))

res1 = minimize_scalar(RotCompactness, args=s1_loads, bounds=(-np.pi/4, 0), method='bounded')
a1 = res1.x; R1 = np.eye(nfac); R1[0,:] = [np.cos(a1), -np.sin(a1)]; R1[1,:] = [np.sin(a1), np.cos(a1)]

res2 = minimize_scalar(RotCompactness, args=s2_loads, bounds=(-np.pi/4, 0), method='bounded')
a2 = res2.x; R2 = np.eye(nfac); R2[0,:] = [np.cos(a2), -np.sin(a2)]; R2[1,:] = [np.sin(a2), np.cos(a2)]

s1_factors_rot = s1_factors.deepcopy(); s1_factors_rot.data = np.matmul(s1_factors.data.T,R1).T
s1_loadings_rot = s1_loadings.deepcopy(); s1_loadings_rot.data = np.matmul(s1_loadings.data.T,R1.T).T

s2_factors_rot = s2_factors.deepcopy(); s2_factors_rot.data = np.matmul(s1_factors.data.T,R1).T
s2_loadings_rot = s2_loadings.deepcopy(); s2_loadings_rot.data = np.matmul(s2_loadings.data.T,R2.T).T

allFacs_rot = hs.signals.Signal1D(np.random.random((4, 2048)),signal_type="EDS_TEM")
allFacs_rot.inav[0:nfac] = s1_factors_rot; allFacs_rot.inav[nfac:nfac+2] = s2_factors_rot
allFacs_rot = setCalibration(allFacs_rot, cal)

hs.plot.plot_spectra(allFacs_rot,style='cascade')
plt.title('Rotated factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_images([s1_loadings_rot, s2_loadings_rot],suptitle='Rotated Loadings', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[1], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                              'left': 0.05, 'right':0.78})


