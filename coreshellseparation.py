# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:31:09 2021

@author: Martin
"""
import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt

core_spect = hs.load('./Spectra/core.msa').data
shell_spect = hs.load('./Spectra/shell.msa').data

size = 50

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

plt.imshow(l[:,:,1])

data_true = l@c
data = np.random.poisson(data_true)
s=hs.signals.Signal1D(data)

s.set_signal_type("EDS_TEM")
s.axes_manager[-1].name = 'E'
s.axes_manager['E'].units = 'keV'
s.axes_manager['E'].scale = 0.010000 #0.01 original
s.axes_manager['E'].offset = -0.2 #-0.02 original
s.set_microscope_parameters(beam_energy=300)

s.change_dtype('float64')

#%% PCA analysis

s.decomposition(normalize_poissonian_noise=True)
s.plot_explained_variance_ratio(n=20, vline=True)

nfac = 2 #antalet faktorer att ta med

#Residualer, kan skippas
sc1 = s.get_decomposition_model(nfac)
res1 = (s - sc1)

res_map=res1.sum(-1).data
res_spec=res1.sum([0,1]).data

#Välj nfac antal komponenter att gå vidare med
s1_factors = s.get_decomposition_factors().inav[0:nfac]
s1_loadings = s.get_decomposition_loadings().inav[0:nfac]

plt.imshow(s1_loadings.data[1,:,:])
plt.plot(s1_factors.data[1,:])

#%% Manuell test av rotation
# Förutsätter att kärnan dominerar och att PCA komponent 2 är skalet
# Value ger ett värde på hur "kompakt" kärnan blir efter rotation
theta = np.linspace(0, -np.pi/4, 100)

Rotated_core = np.empty((50,50))
value = []

for angle in theta:
    temp = s1_loadings.data[0,:,:]*np.cos(angle) + s1_loadings.data[1,:,:]*np.sin(angle)
    Rotated_core = np.dstack((Rotated_core, temp))
    value = np.append(value, np.abs(temp).sum().sum())

#%% Riktig minimering av kärnans "kompakthet"
# Vinkeln som optimeraren hittar används sedan för att rotera loadings och factors
# Begränsningarna är fortfarande att PC1 är summan och PC2 är differensen för skalet
# samt att rotationen är ganska begränsad, dvs. hyfsat tunnt skal.
# Roteringen kunde gjorts mkt bättre med rotationsmatriser etc..

from scipy.optimize import minimize_scalar
loads = (s1_loadings.data[0,:,:], s1_loadings.data[1,:,:])
facs = (s1_factors.data[0,:], s1_factors.data[1,:])

def RotCompactness(a, *loads):

    return np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1]).sum(axis=(0,1))

res = minimize_scalar(RotCompactness, args=loads, bounds=(-np.pi/4, 0), method='bounded')

a = res.x


#OBS! Rotationen är inte korrekt utförd här, det återstår arbete...
s1_loadings_rot = np.cos(a)*loads[0]+np.sin(a)*loads[1]
s1_loadings_rot = np.dstack((s1_loadings_rot, np.cos(a)*loads[1]-np.sin(a)*loads[0]))

s1_factors_rot = np.cos(a)*facs[0]+np.sin(a)*facs[1]
s1_factors_rot = np.vstack((s1_factors_rot, np.cos(a)*facs[1]-np.sin(a)*facs[0]))

plt.plot(s1_factors_rot[0,:])
plt.imshow(s1_loadings_rot[:,:,0])