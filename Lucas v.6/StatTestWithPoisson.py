# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 11:28:17 2021

@author: Lucas
"""

# %matplotlib qt
import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np
from coreshellp import CoreShellP, CoreShellSpec
#from edxmat import EdxMat

dens = 1
x = CoreShellP(50,20.0,15.0,dens,dens,1) #50x50 pixlar,r = 20nm (d=40),core r = 15nm, "densities=1nm^-3", pixel length 1nm

sAg = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
sCu = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM") 
a = CoreShellSpec(x,sCu,sAg)
particle = a.core + a.shell

p = hs.signals.Signal1D(particle)#Läser som HyperSpy-signal
p.add_poissonian_noise()#keep_dtype=True);#Lägger till Poisson-brus

# Filips Metadata:
p.set_signal_type("EDS_TEM") 
p.axes_manager.signal_axes[0].units = 'keV' #OBS! Enheten på energiaxeln är viktig för att kunna plotta
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.add_elements(['Ag','Cu']) #Lägger in ämnen igen tydligen förs de inte med 
p.add_lines(['Ag_La','Cu_Ka'])

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
p.get_calibration_from(cal)
p.axes_manager.indices = [25,25] #Sätter "pekaren" på mitten av partikeln, för att få spektrat från den punkten i plotten
#p.plot(True)  

im = p.get_lines_intensity()

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})

# p.decomposition(True,algorithm="sklearn_pca") # Samma som s.decomposition(normalize_poissonian_noise=True,algorithm="sklearn_pca")
# # #Scree plot för att avgöra rimlig mängd faktorer
# p.plot_explained_variance_ratio(n=50)


# %% Baserat på Nicoletti et al. 
'''
(...) the spectrum images were scaled to normalize the Poisson noise [according to Keenan] and factorized using a 
projected gradient method NMF algorithm. The NMF was repeated for different numbers of components, ranging from four 
to twelve, and showed that eight components were optimal.  
'''


for i in range(6,7):
    dim = i
    #p.decomposition(True,algorithm="sklearn_pca",output_dimension =dim)
    p.decomposition(True,algorithm='NMF',output_dimension =dim)
    #p.plot_decomposition_results()
    factors = p.get_decomposition_factors() #Tar ut faktorerna dvs spektrum
    
    for f in factors:

        f.data /= f.data.max()
    
    loadings =p.get_decomposition_loadings() #loadings är det återskapade bilderna baserat på faktorerna 
    hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
    plt.title('NMF'+str(i))
    plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
    plt.axvline(8040, c='k', ls=':', lw=0.5)
    plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
    plt.axvline(930, c='k', ls=':', lw=0.5)
    plt.axvline(2984, c='k', ls=':', lw=0.5)
    plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    
    hs.plot.plot_images(loadings, cmap='mpl_colors',
                        axes_decor='off', per_row=3,
                        scalebar=[0], scalebar_color='white',
                        padding={'top': 0.95, 'bottom': 0.05,
                                  'left': 0.05, 'right':0.78})
    
    p.blind_source_separation(number_of_components=dim)#,algorithm="orthomax"
    bssfac = p.get_bss_factors()
    bssload = p.get_bss_loadings()
    
    hs.plot.plot_spectra(bssfac.isig[0.0:10000.0],style='cascade') 
    plt.title('NMF+BSS with fastICA')
    plt.axvline(8040, c='k', ls=':', lw=0.5)
    plt.text(x=8040, y=1.6, s='Cu-K$_\\alpha$', color='k')
    plt.axvline(2984, c='k', ls=':', lw=0.5)
    plt.text(x=2984, y=1.6, s='Ag-L$_\\alpha$', color='k')
    plt.axvline(930, c='k', ls=':', lw=0.5)
    plt.text(x=930, y=1.6, s='Cu-L$_\\alpha$', color='k')


    hs.plot.plot_images(bssload, cmap='mpl_colors',
                        axes_decor='off', per_row=3,
                scalebar=[0], scalebar_color='white',
                padding={'top': 0.95, 'bottom': 0.05,
                      'left': 0.05, 'right':0.78})




# %%


