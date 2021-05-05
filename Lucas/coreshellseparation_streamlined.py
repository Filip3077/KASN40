# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 19:21:47 2021

@author: Martin
"""
import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt


#Metod för att skapa STEM-EDX data-kub
core_spect = hs.load('../Spectra/core.msa').data
shell_spect = hs.load('../Spectra/shell.msa').data

size = 50

ld = np.abs(np.linspace(0,size-1,size).reshape((size,1))@np.ones(size).reshape((1,size))-(size-1)/2)
ld = np.sqrt(ld**2 + ld.T**2)

d_core = np.ones((50,50))*15
d_shell = 5

l_core = np.sqrt(d_core**2 - ld**2)
l_core = np.nan_to_num(l_core)

l_shell = np.sqrt((d_core+d_shell )**2- ld**2)-l_core
l_shell = np.nan_to_num(l_shell)

c = np.ones((size,2,2048))
c[:,0,:] = core_spect
c[:,1,:] = shell_spect

l=np.transpose(np.array([l_core, l_shell]))*0.03 # signalstyrkan styrs här

#plt.imshow(l[:,:,1])

data_true = l@c
data = np.random.poisson(data_true)
s=hs.signals.Signal1D(data)

#Mått på signalstyrka
sumsignal = np.sum(data, axis=(0,1,2))
meansignal = np.sum(data, axis=(0,1,2))/size**2

print("Total counts: "+str(sumsignal))
print("Mean counts: "+str(meansignal))

s.set_signal_type("EDS_TEM")
s.axes_manager[-1].name = 'E'
s.axes_manager['E'].units = 'keV'
s.axes_manager['E'].scale = 0.010000 #0.01 original
s.axes_manager['E'].offset = -0.2 #-0.02 original
s.set_microscope_parameters(beam_energy=300)

#Binning och omvandla datatyp
newbin = [1, 1, 8]
s.change_dtype('float64')
s = s.rebin(scale = newbin)

#Energi-skala för plottning
keV = np.linspace(-0.2, -0.2+0.01*2048, int(2048/newbin[2]))
keV_orig = np.linspace(-0.2, -0.2+0.01*2048, 2048)

#%% PCA analysis

s.decomposition(normalize_poissonian_noise=True)
nfac = 2 #antalet faktorer att ta med

#Residualer, kan skippas
sc1 = s.get_decomposition_model(nfac)
res1 = (s - sc1)

res_map=res1.sum(-1).data
res_spec=res1.sum([0,1]).data

#Välj nfac antal komponenter att gå vidare med
s1_factors = [s.get_decomposition_factors().inav[0].data, s.get_decomposition_factors().inav[1].data]
s1_loadings = [s.get_decomposition_loadings().inav[0].data, s.get_decomposition_loadings().inav[1].data]

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
ax1.imshow(s1_loadings[0])
ax2.plot(keV, s1_factors[0])
ax2.set_xlim([-0.1, 10])

#%% Redistribution of intensity between loadings

from scipy.optimize import minimize_scalar

def checkshell(loadings):
    #Checkshell ser bara till att 2:a PC har positiv intensitet i skalet
    x, y = loadings[0].shape
    xc = int(x/2)
    yc = int(y/2)
    
    return -np.mean(loadings[1][yc-2:yc+2, xc-2:xc+2])/np.abs(np.mean(loadings[1][yc-2:yc+2, xc-2:xc+2]))

def RotCompactnessMoment(a, *loads):
    #Gör en rotation och ger "vridmomentet" för kärnarn. För en "kompakt"
    #kärna ska detta vara litet eftersom all intensitet är nära mitten.
    return np.sum(np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1])*ld**4, axis=(0,1))

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

def CorrectShell(b, *args):
    temp_shell = args[1]+b*args[0]
    temp_profile = radial_profile(temp_shell, [24.5, 24.5])
    
    return np.sum((temp_profile/np.max(temp_profile)-args[2])**2)

#Test så att principalkomponent 2 har positiv loading i skalet
t = checkshell(s1_loadings)
s1_loadings[1] = s1_loadings[1]*t
s1_factors[1] = s1_factors[1]*t

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
ax1.imshow(s1_loadings[1])
ax2.plot(keV, s1_factors[1])
ax2.set_xlim([-0.1, 10])

#Rotering för att få fram en ren kärn-loading
res = minimize_scalar(RotCompactnessMoment, args=tuple(s1_loadings), bounds=(-np.pi/2, np.pi/2), method='bounded')

a = res.x

s1_loadings_rot = [np.cos(a)*s1_loadings[0]+np.sin(a)*s1_loadings[1], (np.cos(a)*s1_loadings[1]-np.sin(a)*s1_loadings[0]) ]
s1_factors_rot = [np.cos(a)*s1_factors[0]+np.sin(a)*s1_factors[1], (np.cos(a)*s1_factors[1]-np.sin(a)*s1_factors[0]) ]

#Uppskatta storlek på kärna och kärna+skal utifrån deras areor. Enklast möjliga
#metod. Sist skapas en "ideal" radiell profil som visar hur skalets loading
#bör se ut.
dest_tot = np.sqrt(np.sum((s1_loadings[0] >0.01), axis=(0,1))/np.pi)
dest_core = np.sqrt(np.sum((s1_loadings_rot[0] >0.01), axis=(0,1))/np.pi)
dest_shell = dest_tot-dest_core

center_int_frac = dest_shell/np.sqrt(dest_tot**2-dest_core**2)

shell_profile = radial_profile(s1_loadings_rot[1], [24.5, 24.5])

shell_profile_ideal = np.arange(0,shell_profile.shape[0],1)**2
shell_profile_ideal = np.nan_to_num(np.sqrt(dest_tot**2-shell_profile_ideal)) - np.nan_to_num(np.sqrt(dest_core**2-shell_profile_ideal))
shell_profile_ideal= shell_profile_ideal / np.max(shell_profile_ideal)

#Skalet har vid det här laget oftast antingen ett "hål" i mitten, eller är 
#tjockare uppe på partikeln. Kärnans loading läggs till eller dras ifrån skalet
#för att skapa en skal-loading som har "rätt" profil.
res1 = minimize_scalar(CorrectShell, args=tuple(s1_loadings_rot)+(shell_profile_ideal,), bounds=(-0.8, 0.8), method='bounded')

b = res1.x

#Skalets loading korrigeras för att ge rätt profil. Ändringen kompenseras 
#genom att ändra kärnans factor med motsvarande konstant: ex. om skalet hade 
#"hål" från början måste lite silver i detta fall läggas till kärnan för att 
#skalet ska kunna täcka kärnan jämnt.
s1_loadings_corr = [s1_loadings_rot[0], s1_loadings_rot[1]+b*s1_loadings_rot[0]]
s1_factors_corr = [s1_factors_rot[0]-b*s1_factors_rot[1], s1_factors_rot[1] ]

shell_profile2 = radial_profile(s1_loadings_corr[1], [24.5, 24.5])
core_profile = radial_profile(s1_loadings_corr[0], [24.5, 24.5])

#%% Plots

fig, (ax3, ax4) = plt.subplots(nrows=1, ncols=2)
ax3.imshow(s1_loadings_corr[0])
ax4.imshow(s1_loadings_corr[1])

fig, ax2 = plt.subplots()
ax2.plot(shell_profile/np.max(shell_profile), 'C1', label="Rotation")
ax2.plot(shell_profile_ideal, 'C0', label="True")
ax2.plot(shell_profile2/np.max(shell_profile2), 'C2', label="Corrected")
ax2.plot(core_profile/np.max(core_profile), 'C3', label="Core")
ax2.legend()

fig, ax = plt.subplots()
ax.plot(keV, s1_factors_rot[0], 'C1', label="Rotation")
ax.plot(keV, s1_factors_corr[0], 'C2', label="Corrected")
#ax.plot(keV_orig, core_spect*30, 'C0', label="True")
ax.legend()
ax.set_xlim([-0.1, 10])

#%%


