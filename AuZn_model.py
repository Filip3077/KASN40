# -*- coding: utf-8 -*-
"""
Created on Wed May  5 16:02:42 2021

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
#%% Ladda spektra
sAu=hs.load('./Spectra/20 nm cube Zn0Au100 SDD.msa', signal_type='EDS_TEM');
sZn=hs.load('./Spectra/20 nm cube Zn100Au0 SDD.msa', signal_type='EDS_TEM');
Cback=hs.load('./Spectra/Carbonbackground.msa',signal_type='EDS_TEM')
Cuback= hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block of Cu100*.msa",signal_type = "EDS_TEM")
dim=int(input('Decomposition dimension :'))

background=Cback+Cuback;
#%%Skapa partikel
px=CoreShellP(50,15,10,1,1,1)
p=CoreShellSpec(px,sAu,sZn,signal=True)
p.add_background(background,2)
p=p.getmatr();
p.set_signal_type('EDS_TEM')
p.get_calibration_from(sAu);
p.axes_manager[0].name='x';
p.axes_manager[1].name='y';
p.axes_manager[-1].name='E'
p.axes_manager['E'].unit='eV'
p.add_poissonian_noise()
p.map(gaussian_filter,sigma=2.0)
cut_spectrum_bottom(p,900)
cut_spectrum_range(p,7800,8500)
#%%NMF
p.decomposition(True, algorithm='NMF',output_dimension=dim)
factors=p.get_decomposition_factors();
loadings=p.get_decomposition_loadings();

hs.plot.plot_images(loadings, cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
plt.title('NMF Loadings of simulated Au@Zn particles')
hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
plt.text(x=8640, y=0.8, s='Zn-K$_\\alpha$', color='k')
plt.axvline(8640, c='k', ls=':', lw=0.5)
plt.text(x=1012, y=0.8, s='Zn-L$_\\alpha$', color='k')
plt.axvline(1012, c='k', ls=':', lw=0.5)
plt.axvline(9713, c='k', ls=':', lw=0.5)
plt.text(x=9713, y=0.8, s='Ag-M$_\\alpha$', color='k')
plt.title('NMF factors of simulated Au@Zn particles')


