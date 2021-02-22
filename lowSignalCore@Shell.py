# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:50:02 2021

@author: Filip
"""

import hyperspy.api as hs
from coreshellp import CoreShellP
from addSpectrum import addSpectrum
import numpy as np
from specMapDiff import specMapDiff,cLoadsFacs

s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]
sC = hs.load("./Spectra/Carbonbackground.msa",signal_type="EDS_TEM")
'''


'''
size = 200
cs_mat = CoreShellP(size,30.0,20.0,1,1,1)
cs_mat.scale(1e-9)
csShell = addSpectrum(cs_mat.shell,sAg,0.5e-6)

csCore = addSpectrum(cs_mat.core,sCu,0.5e-6)

carbonMat = np.ones((size,size))
sCarbon = addSpectrum(carbonMat,sC,1)

cs = csCore+csShell + sCarbon
cs.add_poissonian_noise(keep_dtype=True)

cs = cs.rebin(scale = [3,3,1])
cs.axes_manager[0].name = 'y'
cs.axes_manager[1].name = 'x'
cs.axes_manager['x'].units = 'nm'
cs.axes_manager['y'].units = 'nm'
cs.axes_manager[-1].name = 'E'
cs.add_elements(['Ag','Cu','C']) #Lägger in element igen tydligen förs de inte med 
cs.add_lines(['Ag_La','Cu_Ka','C_Ka'])


im = cs.get_lines_intensity()
hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})

decomp_dim = 3
cs.decomposition(output_dimension = decomp_dim ,algorithm='NMF')
NMF_facs = cs.get_decomposition_factors()
NMF_loads = cs.get_decomposition_loadings()

cs.blind_source_separation(number_of_components=decomp_dim)#,algorithm="orthomax"
bss_facs = cs.get_bss_factors()
bss_loads = cs.get_bss_loadings()

#for f in NMF_facs:
#    f.data /= f.data.max()

hs.plot.plot_spectra(NMF_facs.isig[0.0:10000.0],style='cascade') 
hs.plot.plot_images(NMF_loads, cmap='mpl_colors',
            label=['Cu Core', 'Ag Shell'],
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
    
    
hs.plot.plot_spectra(bss_facs.isig[0.0:10000.0],style='cascade') 
hs.plot.plot_images(bss_loads, cmap='mpl_colors',
            axes_decor='off', per_row=1,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})

    
avgcounts = cs.inav[:,:].data.sum()/(cs.data.shape[0]*cs.data.shape[1])
print("\n")
print("Medelantal counts: "+str(avgcounts))

kfacs = [1,0.72980399]
bss_facs.add_elements(['Ag','Cu'])
bss_facs.add_lines(['Ag_La','Cu_Ka'])
print("\n")
print("BSS Ag wt%")
for spectrum in bss_facs:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    print(spectrum.quantification(intensities, method='CL', factors=kfacs,composition_units='weight')[1].data[0])
    
print("\n")
print("NMF Ag wt% (OBS! annorlunda ordning än i BSS)")
NMF_facs.add_elements(['Ag','Cu'])
NMF_facs.add_lines(['Ag_La','Cu_Ka'])
for spectrum in NMF_facs:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    print(spectrum.quantification(intensities, method='CL', factors=kfacs,composition_units='weight')[1].data[0])

NMF = cLoadsFacs(NMF_loads,NMF_facs)
NMF.set_signal_type("EDS_TEM")
NMF.get_calibration_from(cs)


csCore = csCore.rebin(scale = [3,3,1])
csShell = csShell.rebin(scale = [3,3,1])
sCarbon = sCarbon.rebin(scale = [3,3,1])
NMF_bg = specMapDiff(NMF.inav[0],sCarbon)
NMF_1 = specMapDiff(NMF.inav[1],csCore)
NMF_2 =specMapDiff(NMF.inav[2],csShell)

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})

