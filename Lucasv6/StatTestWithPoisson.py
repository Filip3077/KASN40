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

dens = 1
x = CoreShellP(50,20.0,15.0,dens,dens,1) #50x50 pixlar,r = 20nm (d=40),core r = 15nm, "densities=1nm^-3", pixel length 1nm

sAgp = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCup = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM") 
sAg = 0.88*sAgp.data + 0.12*sCup.data
sCu = 0.9*sCup.data + 0.1*sAgp.data
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


for i in range(2,3):
    dim = i
    #p.decomposition(True,algorithm="sklearn_pca",output_dimension =dim)
    p.decomposition(True,algorithm='NMF',output_dimension =dim)
    #p.plot_decomposition_results()
    factors = p.get_decomposition_factors() #Tar ut faktorerna dvs spektrum
    
    # for f in factors:

    #     f.data /= f.data.max()
    
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




# %% Kvantifikation

core_spec = np.sum(a.core,(0,1))
shell_spec = np.sum(a.shell,(0,1))


co = hs.signals.Signal1D(core_spec)
sh = hs.signals.Signal1D(shell_spec)
co.set_signal_type("EDS_TEM")
sh.set_signal_type("EDS_TEM")
co.axes_manager.signal_axes[0].units = 'keV'
sh.axes_manager.signal_axes[0].units = 'keV'
co.add_elements(['Ag','Cu'])
sh.add_elements(['Ag','Cu'])
co.add_lines(['Ag_La','Cu_Ka'])
sh.add_lines(['Ag_La','Cu_Ka'])

cal = hs.load("../Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
co.get_calibration_from(cal)
sh.get_calibration_from(cal)


kfactors = [1.0, 0.72980399]
bw = co.estimate_background_windows(line_width=[5.0, 2.0])
intensities = co.get_lines_intensity(background_windows=bw)
weight_percent_cu_core = co.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
weight_percent_ag_core = co.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

bw = sh.estimate_background_windows(line_width=[5.0, 2.0])
intensities = sh.get_lines_intensity(background_windows=bw)
weight_percent_cu_shell = sh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
weight_percent_ag_shell = sh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

NMFfac = factors.data
NMFco = hs.signals.Signal1D(NMFfac[0])
NMFsh = hs.signals.Signal1D(NMFfac[1])
NMFco.set_signal_type("EDS_TEM")
NMFsh.set_signal_type("EDS_TEM")
NMFco.axes_manager.signal_axes[0].units = 'keV'
NMFsh.axes_manager.signal_axes[0].units = 'keV'
NMFco.add_elements(['Ag','Cu'])
NMFsh.add_elements(['Ag','Cu'])
NMFco.add_lines(['Ag_La','Cu_Ka'])
NMFsh.add_lines(['Ag_La','Cu_Ka'])

NMFco.get_calibration_from(cal)
NMFsh.get_calibration_from(cal)

bw = NMFco.estimate_background_windows(line_width=[5.0, 2.0])
intensities = NMFco.get_lines_intensity(background_windows=bw)
NMFweight_percent_cu_core = NMFco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
NMFweight_percent_ag_core = NMFco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

bw = NMFsh.estimate_background_windows(line_width=[5.0, 2.0])
intensities = NMFsh.get_lines_intensity(background_windows=bw)
NMFweight_percent_cu_shell = NMFsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
NMFweight_percent_ag_shell = NMFsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

BSSfac = bssfac.data
BSSco = hs.signals.Signal1D(BSSfac[0])
BSSsh = hs.signals.Signal1D(BSSfac[1])
BSSco.set_signal_type("EDS_TEM")
BSSsh.set_signal_type("EDS_TEM")
BSSco.axes_manager.signal_axes[0].units = 'keV'
BSSsh.axes_manager.signal_axes[0].units = 'keV'
BSSco.add_elements(['Ag','Cu'])
BSSsh.add_elements(['Ag','Cu'])
BSSco.add_lines(['Ag_La','Cu_Ka'])
BSSsh.add_lines(['Ag_La','Cu_Ka'])

BSSco.get_calibration_from(cal)
BSSsh.get_calibration_from(cal)

bw = BSSco.estimate_background_windows(line_width=[5.0, 2.0])
intensities = NMFco.get_lines_intensity(background_windows=bw)
BSSweight_percent_cu_core = BSSco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
BSSweight_percent_ag_core = BSSco.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

bw = BSSsh.estimate_background_windows(line_width=[5.0, 2.0])
intensities = BSSsh.get_lines_intensity(background_windows=bw)
BSSweight_percent_cu_shell = BSSsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[1].data
BSSweight_percent_ag_shell = BSSsh.quantification(intensities, method='CL', factors=kfactors,composition_units='weight')[0].data

# print(
#       'NMF:'/n
#       'Cu in core: '+str(NMFweight_percent_cu_core)'%'/n
#       )


# fig, axarr = plt.subplots(3,1)
# s1 = relBSScompTot
# s2 = relBSScompCore
# s3 = relBSScompShell

# hs.plot.plot_images(s1, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
#     colorbar='single', vmin='1th', vmax='99th', scalebar='all',
#     scalebar_color='black', suptitle_fontsize=16,
#     padding={'top':0.8, 'bottom':0.10, 'left':0.05,
#               'right':0.85, 'wspace':0.20, 'hspace':0.10},
#     fig=fig)
# hs.plot.plot_images(s2, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
#     colorbar='single', vmin='1th', vmax='99th', scalebar='all',
#     scalebar_color='black', suptitle_fontsize=16,
#     padding={'top':0.8, 'bottom':0.10, 'left':0.05,
#               'right':0.85, 'wspace':0.20, 'hspace':0.10},
#     fig=fig)
# hs.plot.plot_images(s3, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
#     colorbar='single', vmin='1th', vmax='99th', scalebar='all',
#     scalebar_color='black', suptitle_fontsize=16,
#     padding={'top':0.8, 'bottom':0.10, 'left':0.05,
#               'right':0.85, 'wspace':0.20, 'hspace':0.10},
#     fig=fig)
# fig.canvas.draw()

# fig, axarr = plt.subplots(1,2)

# s1 = hs.signals.Signal1D(scipy.misc.ascent()[100:160:10])

# s2 = hs.signals.Signal1D(scipy.misc.ascent()[200:260:10])

# hs.plot.plot_spectra(s1,

#                         style='cascade',

#                         color=[plt.cm.RdBu(i/float(len(s1)-1))

#                                for i in range(len(s1))],

#                         ax=axarr[0],

#                         fig=fig)

# hs.plot.plot_spectra(s2,

#                         style='cascade',

#                         color=[plt.cm.summer(i/float(len(s1)-1))

#                                for i in range(len(s1))],

#                         ax=axarr[1],

#                         fig=fig)

# axarr[0].set_xlabel('RdBu (colormap)')

# axarr[1].set_xlabel('summer (colormap)')

# fig.canvas.draw()

