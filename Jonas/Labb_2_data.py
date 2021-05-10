# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 08:19:53 2021

Data till labb utförd 22/4 2021 - EDX på 22 och 30 nm Au@Zn-partiklar.

@author: Jonas
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
from core_shell_sim.sim.coreshellp import *
from core_shell_sim.proc.specerr import *
from core_shell_sim.proc.specMapDiff import *
import numpy as np
from core_shell_sim.proc.loadassign import checkLoadFit
from scipy.ndimage import gaussian_filter
from core_shell_sim.qual.k_factors import silver_k_factor
from core_shell_sim.rot.radial_profile import transfer_elements


p1=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\NiCoO\\NiCoOx site 1.rpl',signal_type='EDS_TEM').T
p2=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\NiCoO\\NiCoOx site 2.rpl',signal_type='EDS_TEM').T
p3=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\NiCoO\\NiCoOx site 3.rpl',signal_type='EDS_TEM').T
p4=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\NiCoO\\NiCoOx site 4.rpl',signal_type='EDS_TEM').T
p5=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\NiCoO\\NiCoOx site 5.rpl',signal_type='EDS_TEM').T
ps=[p1,p2,p3,p4,p5];
#%%Skära i rymddimensionen
ps[0]=ps[0].inav[50:230,50:200];
ps[1]=ps[1].inav[55:200,25:175];
ps[2]=ps[2].inav[125:256,47:210]
ps[3]=ps[3].inav[75:200,30:170];
ps[4]= ps[4].inav[70:200,65:210];
#%%Behövliga listor
ims=[];
factors=[];
loadings=[];
#%%Axlar, grundämne
for i in range(len(ps)):
    ps[i].axes_manager[0].name = 'y'
    ps[i].axes_manager[1].name = 'x'
    ps[i].axes_manager['x'].units = 'nm'
    ps[i].axes_manager['y'].units = 'nm'
    ps[i].axes_manager[-1].name = 'E'
    ps[i].axes_manager['E'].units='eV'
    ps[i].axes_manager['E'].offset=-200;#Kan komma att ändras
    ps[i].axes_manager['E'].scale=10;
    #ps[i].axes_manager['x'].scale=0.1666;
    #ps[i].axes_manager['y'].scale=0.1666;
    ps[i].add_elements(['Ni','Co','O','C','Cu'])
    ps[i].add_lines(['O_Ka','Co_Ka','Ni_Ka','C_Ka','Cu_Ka'])
    ps[i].change_dtype('float64')
    ps[i]=ps[i].rebin(scale=[2,2,1])
    ims.append(ps[i].get_lines_intensity())
for p in ps:
    p=p.map(gaussian_filter,sigma=2.0)

    
#%%PCA 1
ps[0].decomposition('sklearn_pca')
ps[0].plot_explained_variance_ratio();
#%%NMF 1
ps[0].decomposition(True, algorithm='NMF', output_dimension=4)
if len(factors)>0:
    factors[0]=ps[0].get_decomposition_factors()
else:
    factors.append(ps[0].get_decomposition_factors());
if len(loadings)>0:
    loadings[0]=ps[0].get_decomposition_loadings()
else:
    loadings.append(ps[0].get_decomposition_loadings())
#%%Plotta 1

redBlueMap(ims[0],'Site 1')
for i in range(len(factors[0])):
    factors[0].inav[i].data/=factors[0].inav[i].data.max()
hs.plot.plot_spectra(factors[0].isig[0.0:10000.0],style='cascade')
plt.title('Site 1')

hs.plot.plot_images(loadings[0],suptitle='Site 1', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 2
ps[1].decomposition('sklearn_pca')
ps[1].plot_explained_variance_ratio();
#%%NMF 2
ps[1].decomposition(True, algorithm='NMF', output_dimension=4)
if len(factors)>1:
    factors[1]=ps[1].get_decomposition_factors()
else:
    factors.append(ps[1].get_decomposition_factors());
if len(loadings)>1:
    loadings[1]=ps[1].get_decomposition_loadings()
else:
    loadings.append(ps[1].get_decomposition_loadings())
#%%Plotta 2

redBlueMap(ims[1],'Site 2')
for i in range(len(factors[1])):
    factors[1].inav[i].data/=factors[1].inav[i].data.max()
hs.plot.plot_spectra(factors[1].isig[0.0:10000.0],style='cascade')
plt.title('Site 2')

hs.plot.plot_images(loadings[1],suptitle='Site 2', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 3
ps[2].decomposition('sklearn_pca')
ps[2].plot_explained_variance_ratio();
#%%NMF 3
ps[2].decomposition(True, algorithm='NMF', output_dimension=4)
if len(factors)>2:
    factors[2]=ps[2].get_decomposition_factors()
else:
    factors.append(ps[2].get_decomposition_factors());
if len(loadings)>2:
    loadings[2]=ps[2].get_decomposition_loadings()
else:
    loadings.append(ps[2].get_decomposition_loadings())
#%%Plotta 3

redBlueMap(ims[2],'Site 2')
for i in range(len(factors[2])):
    factors[2].inav[i].data/=factors[2].inav[i].data.max()

hs.plot.plot_spectra(factors[2].isig[0.0:10000.0],style='cascade')
plt.title('Site 2')

hs.plot.plot_images(loadings[2],suptitle='Site 3', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%PCA 4
ps[3].decomposition('sklearn_pca')
ps[3].plot_explained_variance_ratio();
#%%NMF 4
ps[3].decomposition(True, algorithm='NMF', output_dimension=4)
if len(factors)>3:
    factors[3]=ps[3].get_decomposition_factors()
else:
    factors.append(ps[3].get_decomposition_factors());
if len(loadings)>3:
    loadings[3]=ps[3].get_decomposition_loadings()
else:
    loadings.append(ps[3].get_decomposition_loadings())
#%%Plotta 4

redBlueMap(ims[3],'Site 4')
for i in range(len(factors[3])):
    factors[3].inav[i].data/=factors[3].inav[i].data.max()
hs.plot.plot_spectra(factors[3].isig[0.0:10000.0],style='cascade')
plt.title('Site 4')

hs.plot.plot_images(loadings[3],suptitle='Site 4', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 5
ps[4].decomposition('sklearn_pca')
ps[4].plot_explained_variance_ratio();
#%%NMF 5
ps[4].decomposition(True, algorithm='NMF', output_dimension=2)
if len(factors)>4:
    factors[4]=ps[4].get_decomposition_factors()
else:
    factors.append(ps[4].get_decomposition_factors());
if len(loadings)>4:
    loadings[4]=ps[4].get_decomposition_loadings()
else:
    loadings.append(ps[4].get_decomposition_loadings())
#%%Plotta 5

redBlueMap(ims[4],'Site 5')
for i in range(len(factors[4])):
    factors[4].inav[i].data/=factors[4].inav[i].data.max()

hs.plot.plot_spectra(factors[4].isig[0.0:10000.0],style='cascade')
plt.title('Site 5')

hs.plot.plot_images(loadings[4],suptitle='Site 5', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%Kvantifiera sista partikeln
quant=[]
for f in factors:
    f.set_signal_type("EDS_TEM")
    #f.axes_manager.signal_axes[0].units='eV'
    f.get_calibration_from(ps[4])
    f.add_elements(['Ni','Co','O'])#,'Cu','C'
    f.add_lines(['Ni_Ka','Co_Ka','O_Ka'])#,'Cu_Ka','C_Ka'
    #bg = f.estimate_background_windows(line_width=[3.0, 5.0])
    intensities=f.get_lines_intensity()
    quant.append(f.quantification(intensities,method='CL',factors=[1.222,1.219,1.903],composition_units='weight'))
