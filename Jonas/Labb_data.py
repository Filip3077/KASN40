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

p22_1=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 22 nm\\EDS Data_particle 1.rpl',signal_type='EDS_TEM').T
p22_2=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 22 nm\\EDS Data_particle 2.rpl',signal_type='EDS_TEM').T
p22_3=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 22 nm\\EDS Data_particle 3.rpl',signal_type='EDS_TEM').T
p22=[p22_1,p22_2,p22_3];
p30_1=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 1.rpl',signal_type='EDS_TEM').T
p30_2=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 2.rpl',signal_type='EDS_TEM').T
p30_3=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 3.rpl',signal_type='EDS_TEM').T
p30_4=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 4.rpl',signal_type='EDS_TEM').T
p30=[p30_1,p30_2,p30_3,p30_4];
#%%Behövliga listor
im22=[];
im30=[];
factors22=[];
loadings22=[];
factors30=[];
loadings30=[];
#%%Axlar, grundämne
for i in range(len(p22)):
    p22[i].axes_manager[0].name = 'y'
    p22[i].axes_manager[1].name = 'x'
    p22[i].axes_manager['x'].units = 'nm'
    p22[i].axes_manager['y'].units = 'nm'
    p22[i].axes_manager[-1].name = 'E'
    p22[i].axes_manager['E'].units='eV'
    p22[i].axes_manager['E'].offset=-200;
    p22[i].axes_manager['E'].scale=10;
    p22[i].axes_manager['x'].scale=0.1666;
    p22[i].axes_manager['y'].scale=0.1666;
    #cut_spectrum_bottom(p22[i],900)
    #cut_spectrum_range(p22[i],7900,8400)
    p22[i].add_elements(['Au','Zn','C','Cu'])
    p22[i].add_lines(['Au_Ma','Zn_La','C_Ka','Cu_Ka'])
    p22[i].change_dtype('float64')
    p22[i]=p22[i].rebin(scale=[5,5,1])
    im22.append(p22[i].get_lines_intensity())
for p in p22:
    p=p.map(gaussian_filter,sigma=3.0)
for i in range(len(p30)):
    p30[i].axes_manager[0].name = 'y'
    p30[i].axes_manager[1].name = 'x'
    p30[i].axes_manager['x'].units = 'nm'
    p30[i].axes_manager['y'].units = 'nm'
    p30[i].axes_manager[-1].name = 'E'
    p30[i].axes_manager['E'].units='eV'
    p30[i].axes_manager['E'].scale=10;
    p30[i].axes_manager['E'].offset=170;
    #cut_spectrum_bottom(p30[i],900.0)
    #cut_spectrum_range(p30[i],7900.0,8400.0)
    p30[i].add_elements(['Au','Zn','C','Cu'])
    p30[i].add_lines(['Au_Ma','Zn_La','C_Ka','Cu_Ka'])
    p30[i].change_dtype('float64')
    p30[i]=p30[i].rebin(scale=[5,5,1])
    im30.append(p30[i].get_lines_intensity());
for p in p30:
    p=p.map(gaussian_filter,sigma=3.0)
    
#%%PCA 22_1
p22[0].decomposition('sklearn_pca')
p22[0].plot_explained_variance_ratio();
#%%NMF 22_1
p22[0].decomposition(True, algorithm='NMF', output_dimension=2)
if len(factors22)>0:
    factors[0]=p22[0].get_decomposition_factors()
else:
    factors22.append(p22[0].get_decomposition_factors());
if len(loadings22)>0:
    loadings22[0]=p22[0].get_decomposition_loadings()
else:
    loadings22.append(p22[0].get_decomposition_loadings())
#%%Plotta 22_1

redBlueMap(im22[0],'20 nm particle 1')

hs.plot.plot_spectra(factors22[0].isig[0.0:10000.0],style='cascade')
plt.title('20 nm particle 1')

hs.plot.plot_images(loadings22[0],suptitle='20 nm particle 1', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 22_2
p22[1].decomposition('sklearn_pca')
p22[1].plot_explained_variance_ratio();
#%%NMF 22_2
p22[1].decomposition(True, algorithm='NMF', output_dimension=3)
if len(factors22)>1:
    factors[1]=p22[1].get_decomposition_factors()
else:
    factors22.append(p22[1].get_decomposition_factors());
if len(loadings22)>1:
    loadings22[1]=p22[1].get_decomposition_loadings()
else:
    loadings22.append(p22[1].get_decomposition_loadings())
#%%Plotta 22_2

redBlueMap(im22[1],'20 nm particle 2')

hs.plot.plot_spectra(factors22[1].isig[0.0:10000.0],style='cascade')
plt.title('20 nm particle 2')

hs.plot.plot_images(loadings22[1],suptitle='20 nm particle 2', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 22_3
p22[2].decomposition('sklearn_pca')
p22[2].plot_explained_variance_ratio();
#%%NMF 22_3
p22[2].decomposition(True, algorithm='NMF', output_dimension=3)
if len(factors22)>2:
    factors22[2]=p22[2].get_decomposition_factors()
else:
    factors22.append(p22[2].get_decomposition_factors());
if len(loadings22)>2:
    loadings22[2]=p22[2].get_decomposition_loadings()
else:
    loadings22.append(p22[2].get_decomposition_loadings())
#%%Plotta 22_3

redBlueMap(im22[2],'20 nm particle 2')

hs.plot.plot_spectra(factors22[2].isig[0.0:10000.0],style='cascade')
plt.title('20 nm particle 2')

hs.plot.plot_images(loadings22[2],suptitle='20 nm particle 3', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%PCA 30_1
p30[0].decomposition('sklearn_pca')
p30[0].plot_explained_variance_ratio();
#%%NMF 30_1
p30[0].decomposition(True, algorithm='NMF', output_dimension=3)
if len(factors30)>0:
    factors30[0]=p30[0].get_decomposition_factors()
else:
    factors30.append(p30[0].get_decomposition_factors());
if len(loadings30)>0:
    loadings30[0]=p30[0].get_decomposition_loadings()
else:
    loadings30.append(p30[0].get_decomposition_loadings())
#%%Plotta 30_1

redBlueMap(im30[0],'30 nm particle 1')

hs.plot.plot_spectra(factors30[0].isig[0.0:10000.0],style='cascade')
plt.title('30 nm particle 1')

hs.plot.plot_images(loadings30[0],suptitle='30 nm particle 1', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 30_2
p30[1].decomposition('sklearn_pca')
p30[1].plot_explained_variance_ratio();
#%%NMF 30_2
p30[1].decomposition(True, algorithm='NMF', output_dimension=3)
if len(factors30)>1:
    factors[1]=p30[1].get_decomposition_factors()
else:
    factors30.append(p30[1].get_decomposition_factors());
if len(loadings30)>1:
    loadings30[1]=p30[1].get_decomposition_loadings()
else:
    loadings30.append(p30[1].get_decomposition_loadings())
#%%Plotta 30_2

redBlueMap(im30[1],'30 nm particle 2')

hs.plot.plot_spectra(factors30[1].isig[0.0:10000.0],style='cascade')
plt.title('30 nm particle 1')

hs.plot.plot_images(loadings30[1],suptitle='30 nm particle 1', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%PCA 30_3
p30[2].decomposition('sklearn_pca')
p30[2].plot_explained_variance_ratio();
#%%NMF 30_3
p30[2].decomposition(True, algorithm='NMF', output_dimension=3)
if len(factors30)>2:
    factors[2]=p30[2].get_decomposition_factors()
else:
    factors30.append(p30[2].get_decomposition_factors());
if len(loadings30)>2:
    loadings30[2]=p30[2].get_decomposition_loadings()
else:
    loadings30.append(p30[2].get_decomposition_loadings())
#%%Plotta 30_3

redBlueMap(im30[2],'30 nm particle 3')

hs.plot.plot_spectra(factors30[2].isig[0.0:10000.0],style='cascade')
plt.title('30 nm particle 3')

hs.plot.plot_images(loadings30[2],suptitle='30 nm particle 3', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%PCA 30_4
p30[3].decomposition('sklearn_pca')
p30[3].plot_explained_variance_ratio();
#%%NMF 30_4
p30[3].decomposition(True, algorithm='NMF', output_dimension=5)
if len(factors30)>3:
    factors30[3]=p30[3].get_decomposition_factors()
else:
    factors30.append(p30[3].get_decomposition_factors());
if len(loadings30)>3:
    loadings30[3]=p30[3].get_decomposition_loadings()
else:
    loadings30.append(p30[3].get_decomposition_loadings())
#%%Plotta 30_4

redBlueMap(im30[3],'30 nm particle 4')

hs.plot.plot_spectra(factors30[3].isig[0.0:10000.0],style='cascade')
plt.title('30 nm particle 4')

hs.plot.plot_images(loadings30[3],suptitle='30 nm particle 4', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
#%%Kvantifiera partiklar
quant30=[]
for f in factors30:
    f.set_signal_type("EDS_TEM")
    f.get_calibration_from(p30[3])
    f.add_elements(['Au','Zn'])
    f.add_lines(['Au_Ma','Zn_La'])
    bg = f.estimate_background_windows(line_width=[5.0, 7.0])
    intensities=f.get_lines_intensity(background_windows=bg)
    quant30.append(f.quantification(intensities,method='CL',factors=[1.696,1.659],composition_units='weight'))
