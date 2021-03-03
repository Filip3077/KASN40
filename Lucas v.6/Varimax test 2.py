# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:17:23 2021

@author: Lucas
"""

from varimax import varimax
import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellfunctions import *
from specMapDiff import *
import numpy as np
from numpy import linalg

#%% Particle Generation
sAgPure = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")

sAg = 0.8*sAgPure.data + 0.2*sCuPure.data
sCu = 0.9*sCuPure.data + 0.1*sAgPure.data


a = genfullparticle(50,20,15,sCu,sAg)

p = a.full
p.add_poissonian_noise()

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
p = setCalibration(p, cal)

p.change_dtype('float64')

#%% PCA part

p.decomposition(normalize_poissonian_noise=True)#,algorithm="sklearn_pca") #
p.plot_explained_variance_ratio(n=20, vline=True)

nfac = 2 

sc1 = p.get_decomposition_model(nfac)
res1 = (p - sc1)

res_map=res1.sum(-1).data[0,:]
res_spec=res1.sum([0,1]).data

p_factors = p.get_decomposition_factors()
p_loadings = p.get_decomposition_loadings()


hs.plot.plot_spectra(p_factors.inav[0:nfac].isig[0.0:10000.0],style='cascade')
plt.title('SVD factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_images(p_loadings.inav[0:nfac],suptitle='SVD Loadings', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                              'left': 0.05, 'right':0.78})

#%% Varimax part

p_factors_selected = p_factors.inav[0:nfac]
p_loadings_selected = p_loadings.inav[0:nfac]

#Unfold för att göra er 2+1D matris till 1+1D som varimax vill ha.
p_loadings_unfold = p_loadings_selected.deepcopy()
p_loadings_unfold.unfold()

p_factors_selected = p_factors_selected.data
p_loadings_selected = p_loadings_unfold.data

R =varimax(p_loadings_selected.T)

p_loadings_selected_rot = np.matmul(p_loadings_selected.T, R).T
p_factors_selected_rot = np.matmul(linalg.inv(R), p_factors_selected)

#Ibland blir faktorer och loadings negativa. Då är det bara att manuellt flippa dem.
for i in range(nfac):
    if p_factors_selected_rot[i,:].sum() < 0:
        p_factors_selected_rot[i,:] = -p_factors_selected_rot[i,:]
    if p_loadings_selected_rot[i,:].sum() < 0:
        p_loadings_selected_rot[i,:] = -p_loadings_selected_rot[i,:]

#"Vik ihop" matrisen till ursprungligt format, fast nu roterad
p_loadings_unfold.data = p_loadings_selected_rot
p_loadings_unfold.fold()

#%% 
maxint = np.max(p_factors_selected_rot[0:nfac,:])*1.1

keV = np.arange(-0.2, 20.28, 0.01) 
# fig1, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))
fig1, ax1 = plt.subplots()

ax1.plot(keV, p_factors_selected_rot[0,:], label='C1')
ax1.plot(keV, p_factors_selected_rot[1,:], label='C2')
# ax1.plot(keV, p_factors_selected_rot[2,:], label='C3')
# ax1.plot(keV, p_factors_selected_rot[3,:], label='C4')
#ax1.plot(keV, p_factors_selected_rot[4,:], label='C5')
ax1.axis([-0.2,20.28,-0.5, maxint])
ax1.legend()
plt.title('Rotated factors')

fig2, (ax3, ax4) = plt.subplots(nrows=1, ncols=2, figsize=(20, 4))
ax3.imshow(p_loadings_unfold.inav[0].data, cmap=plt.get_cmap('viridis'), vmin=-2, vmax=7)
ax4.imshow(p_loadings_unfold.inav[1].data, cmap=plt.get_cmap('viridis'), vmin=-2, vmax=7)
# ax5.imshow(p_loadings_unfold.inav[2].data, cmap=plt.get_cmap('viridis'), vmin=-2, vmax=7)

#%% Our regular plots

p_fac_rot = hs.signals.Signal1D(p_factors_selected_rot)
p_fac_rot = setCalibration(p_fac_rot, cal)
p_load_rot = hs.signals.BaseSignal(p_loadings_unfold)
ploadrot = p_load_rot.transpose(signal_axes=2)

hs.plot.plot_spectra(p_fac_rot.isig[0.0:10000.0],style='cascade')
plt.title('Varimax Factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_images(ploadrot,suptitle='Varimax Loadings', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                              'left': 0.05, 'right':0.78})
