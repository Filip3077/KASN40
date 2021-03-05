# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 13:33:08 2021

@author: Lucas
"""

from varimax import varimax
import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellfunctions import *
from specMapDiff import *
import numpy as np
from numpy import linalg

#%%
sAgPure = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")

sAg = 0.88*sAgPure.data + 0.12*sCuPure.data
sCu = 0.9*sCuPure.data + 0.1*sAgPure.data


a = genfullparticle(50,20,15,sCu,sAg)

p = a.full
p.add_poissonian_noise()

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
p = setCalibration(p, cal)

p.change_dtype('float64')
#%% PCA part

p.decomposition(normalize_poissonian_noise=True,algorithm="sklearn_pca")
p.plot_explained_variance_ratio(n=20, vline=True)


sc1 = p.get_decomposition_model(2)
res1 = (p - sc1)

res_map=res1.sum(-1).data[0,:]
res_spec=res1.sum([0,1]).data

p_factors = p.get_decomposition_factors()
p_loadings = p.get_decomposition_loadings()


hs.plot.plot_spectra(p_factors.inav[0:2].isig[0.0:10000.0],style='cascade')
plt.title('SVD factors')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

hs.plot.plot_images(p_loadings.inav[0:2],suptitle='SVD Loadings', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                              'left': 0.05, 'right':0.78})


# p_factors.inav[0].isig[0:180].plot() # plot low energy region of factor

#%% Varimax part
sum_loadings = p_loadings.sum(-1)
# p_loadings = sum_loadings
nfac = 2

p_factors_selected = p_factors.data[0:nfac,:]
p_loadings_selected = p_loadings.data[0:nfac,:]
sum_loadings_selected = sum_loadings.data[0:nfac,:]

R = varimax(sum_loadings_selected.T)

p_loadings_selected_rot = np.matmul(sum_loadings_selected.T, R).T
p_factors_selected_rot = np.matmul(linalg.inv(R), p_factors_selected)

p_factors_selected_rot[0,:] = -p_factors_selected_rot[0,:]
p_loadings_selected_rot[0,:] = -p_loadings_selected_rot[0,:]

p_factors_selected_rot[1,:] = -p_factors_selected_rot[1,:]
p_loadings_selected_rot[1,:] = -p_loadings_selected_rot[1,:]

#%% Plotting Something

keV = np.arange(-0.2, 20.28, 0.01) 
# fig1, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))
fig1, ax1 = plt.subplots()

ax1.plot(keV, p_factors_selected_rot[0,:], label='C1')
ax1.plot(keV, p_factors_selected_rot[1,:], label='C2')
# ax1.plot(keV, p_factors_selected_rot[2,:], label='C3')
# ax1.plot(keV, p_factors_selected_rot[3,:], label='C4')
#ax1.plot(keV, p_factors_selected_rot[4,:], label='C5')
ax1.axis([-0.2,20.28,-25, 20])
ax1.legend()
plt.title('Rotated factors')

# ax2.plot(keV, p_factors_selected_rot[0,:], label='C1')
# ax2.plot(keV, p_factors_selected_rot[1,:], label='C2')
# # ax2.plot(keV, p_factors_selected_rot[2,:], label='C3')
# # ax2.plot(keV, p_factors_selected_rot[3,:], label='C4')
# #ax2.plot(keV, p_factors_selected_rot[4,:], label='C5')
# ax2.axis([1.6,12.5,0, 40])
# ax2.legend()

fig2, ax3 = plt.subplots()

ax3.plot(p_loadings_selected_rot[0,:], label='C1')
ax3.plot(p_loadings_selected_rot[1,:], label='C2')
# ax3.plot(p_loadings_selected_rot[2,:], label='C3')
# ax3.plot(p_loadings_selected_rot[3,:], label='C4')
#ax3.plot(p_loadings_selected_rot[4,:], label='C5')
ax3.axis([0,50,-20,20])
ax3.legend()
plt.title('Rotated Loadings')

#%% To our methods

# im = p.get_lines_intensity() # 

# hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
#     colorbar='single', vmin='1th', vmax='99th', scalebar='all',
#     scalebar_color='black', suptitle_fontsize=16,
#     padding={'top':0.8, 'bottom':0.10, 'left':0.05,
#               'right':0.85, 'wspace':0.20, 'hspace':0.10})

# p_factors_selected_rot = setCalibration(hs.signals.Signal1D(p_factors_selected_rot), cal)
# p_loadings_selected_rot = setCalibration(hs.signals.Signal1D(p_loadings_selected_rot), cal)

p_fac_rot = hs.signals.Signal1D(p_factors_selected_rot)
p_fac_rot = setCalibration(p_fac_rot, cal)
p_load_rot = hs.signals.BaseSignal(p_loadings_selected_rot)
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

#%% Keeping the Loadings as pictures

R = varimax(sum_loadings_selected.T)

p_loadings_selected_rot = np.matmul(p_loadings_selected.T, R).T
p_factors_selected_rot = np.matmul(linalg.inv(R), p_factors_selected)

# p_factors_selected_rot[0,:] = -p_factors_selected_rot[0,:]
# p_loadings_selected_rot[0,:] = -p_loadings_selected_rot[0,:]

p_factors_selected_rot[1,:] = -p_factors_selected_rot[1,:]
p_loadings_selected_rot[1,:] = -p_loadings_selected_rot[1,:]

#%% Plotting

keV = np.arange(-0.2, 20.28, 0.01) 
# fig1, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))
fig1, ax1 = plt.subplots()

ax1.plot(keV, p_factors_selected_rot[0,:], label='C1')
ax1.plot(keV, p_factors_selected_rot[1,:], label='C2')
# ax1.plot(keV, p_factors_selected_rot[2,:], label='C3')
# ax1.plot(keV, p_factors_selected_rot[3,:], label='C4')
#ax1.plot(keV, p_factors_selected_rot[4,:], label='C5')
ax1.axis([-0.2,20.28,-10, 20])
ax1.legend()
plt.title('Rotated Factors')

# ax2.plot(keV, p_factors_selected_rot[0,:], label='C1')
# ax2.plot(keV, p_factors_selected_rot[1,:], label='C2')
# # ax2.plot(keV, p_factors_selected_rot[2,:], label='C3')
# # ax2.plot(keV, p_factors_selected_rot[3,:], label='C4')
# #ax2.plot(keV, p_factors_selected_rot[4,:], label='C5')
# ax2.axis([1.6,12.5,0, 40])
# ax2.legend()

fig2, ax3 = plt.subplots()

ax3.plot(p_loadings_selected_rot[0,:], label='C1')
ax3.plot(p_loadings_selected_rot[1,:], label='C2')
# ax3.plot(p_loadings_selected_rot[2,:], label='C3')
# ax3.plot(p_loadings_selected_rot[3,:], label='C4')
#ax3.plot(p_loadings_selected_rot[4,:], label='C5')
ax3.axis([0,50,-2,3])
# ax3.legend()
plt.title('Rotated Loadings')

#%% And usual plotting

p_fac_rot = hs.signals.Signal1D(p_factors_selected_rot)
p_fac_rot = setCalibration(p_fac_rot, cal)
p_load_rot = hs.signals.BaseSignal(p_loadings_selected_rot)
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