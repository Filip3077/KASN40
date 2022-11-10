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

t=int(input('How many reps?:'))
ga=bool(input('Gaussian (1/0)?: '))
p4=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\NiCoO\\NiCoOx site 4.rpl',signal_type='EDS_TEM').T

ps=[p4];
#%%Skära i rymddimensionen
ps[0]=ps[0].inav[75:200,30:170];

#%%Behövliga listor
ims=[];
factors=[];
loadings=[];
quant=[];
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
if ga:
    for p in ps:
        p=p.map(gaussian_filter,sigma=2.0)

    
#%%PCA 4
ps[0].decomposition('sklearn_pca')
ps[0].plot_explained_variance_ratio();
for u in range(t):
#%%NMF 4
    ps[0].decomposition(True, algorithm='NMF', output_dimension=4)

    factors=ps[0].get_decomposition_factors();
    loadings=ps[0].get_decomposition_loadings()

#%%Plotta 4

    
    for i in range(len(factors)):
        factors.inav[i].data/=factors.inav[i].data.max()
    if t==1:
        redBlueMap(ims[0],'Site 4')
        hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
        plt.title('Site 4')
    
        hs.plot.plot_images(loadings,suptitle='Site 4', cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%Kvantifiera sista partikeln

    factors.set_signal_type("EDS_TEM")
    #f.axes_manager.signal_axes[0].units='eV'
    factors.get_calibration_from(ps[0])
    factors.add_elements(['Ni','Co','O'])#,'Cu','C'
    factors.add_lines(['Ni_Ka','Co_Ka','O_Ka'])#,'Cu_Ka','C_Ka'
    #bg = f.estimate_background_windows(line_width=[3.0, 5.0])
    intensities=factors.get_lines_intensity()
    q=factors.quantification(intensities,method='CL',factors=[1.222,1.219,1.903],composition_units='weight')
    quant.append(q)
#%%Skapa statistik10
l2=np.zeros((3,len(quant)));
l4=np.zeros((3,len(quant)));
for i in range(len(quant)):
    for k in range(len(l2)):#Co,Ni,O
        l2[k][i]=quant[i][k].data[1];
        l4[k][i]=quant[i][k].data[3];
mean_2=np.sum(1/t*l2,axis=1)
mean_4=np.sum(1/t*l4,axis=1)
std_2=np.std(l2,axis=1)
std_4=np.std(l4,axis=1)
    