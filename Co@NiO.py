# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 18:07:18 2021

@author: Filip
"""

import hyperspy.api as hs
from core_shell_sim.sim.coreshellp import *
from scipy.ndimage import gaussian_filter

Co = hs.load("./Spectra/Co.msa",signal_type="EDS_TEM")
NiO = hs.load("./Spectra/NiO.msa",signal_type="EDS_TEM")
C = hs.load("./Spectra/Carbonbackground.msa",signal_type = "EDS_TEM")
Cu = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block of Cu100*.msa",signal_type = "EDS_TEM")

background = C + Cu/10


shell = NiO*0.8 + Co*0.2
cs_mat = CoreShellP(170,54,35,1e-2,1e-2,1)
cs = CoreShellSpec(cs_mat,Co,shell)
cs.add_background(background,1)
s = cs.getmatr()
s.add_poissonian_noise(keep_dtype=True)

#s_gauss = gaussian_filter(s.data, sigma=2)
#s = hs.signals.Signal1D(s_gauss)
s_avgcounts = s.inav[:,:].data.sum()/(s.data.shape[0]*s.data.shape[1])
s_totcounts = s.inav[:,:].data.sum()

print('\n')
print("Medelantal counts: "+str(s_avgcounts))
print("Totalt antal counts: " + str(s_totcounts))
print('\n')

#%%
s.decomposition()
s.plot_explained_variance_ratio()
#dim = s.estimate_number_of_clusters(cluster_source="decomposition", metric="gap")
dim = 3
#%%
s.set_signal_type("EDS_TEM") 
s.get_calibration_from(Co)
s.set_elements(['Co','Ni', 'O','Cu','C'])
s.set_lines(['Co_Ka','Ni_Ka','O_Ka','Cu_Ka','C_Ka'])

im1 = s.get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im1,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%%

s.decomposition(algorithm = 'NMF', output_dimension = dim)
s_loads = s.get_decomposition_loadings()
s_facs = s.get_decomposition_factors()

hs.plot.plot_spectra(s_facs,style='cascade') 
hs.plot.plot_images(s_loads, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})




#%%
#BSS

s.blind_source_separation(number_of_components=dim)#,algorithm="orthomax"
bssfac = s.get_bss_factors()
bssload = s.get_bss_loadings()

hs.plot.plot_spectra(bssfac,style='cascade') 
hs.plot.plot_images(bssload, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})

#%%
Ni = hs.load("./Spectra/Ni_true.msa",signal_type="EDS_TEM")

sCal  = Co*0.5 + Ni*0.5
sCal.set_elements(['Co','Ni'])
sCal.set_lines(['Co_Ka','Ni_Ka'])
kfac = [0.8806,1]
bw = sCal.estimate_background_windows(line_width=[1.7, 2.0])
sCal.plot(background_windows=bw)
intensities = sCal.get_lines_intensity(background_windows=bw)

kfac_NiCo = [1/(intensities[0].data[0]/intensities[1].data[0]),1]
atomic_percent = sCal.quantification(intensities, method='CL',
                                  factors=kfac)
print(atomic_percent[1].data[0])


#%%
s_facs.set_signal_type("EDS_TEM") 
s_facs.get_calibration_from(Co)

s_facs.set_elements(['Co','Ni'])
s_facs.set_lines(['Co_Ka','Ni_Ka'])
for each in s_facs.split():
    bw = each.estimate_background_windows(line_width=[1.7, 2.0])
    intensities = each.get_lines_intensity(background_windows=bw)
    atomic_percent = each.quantification(intensities, method='CL',
                                  factors=kfac_NiCo)
    for elem in atomic_percent:
        print(elem.metadata.Sample.elements[0] + ": "+ str(elem.data[0]))
    print("\n")


#%%
bssfac.set_signal_type("EDS_TEM") 
bssfac.get_calibration_from(Co)

bssfac.set_elements(['Co','Ni'])
bssfac.set_lines(['Co_Ka','Ni_Ka'])
for each in bssfac.split():
    bw = each.estimate_background_windows(line_width=[1.7, 2.0])
    intensities = each.get_lines_intensity(background_windows=bw)
    atomic_percent = each.quantification(intensities, method='CL',
                                  factors=kfac_NiCo)
    for elem in atomic_percent:
        print(elem.metadata.Sample.elements[0] + ": "+ str(elem.data[0]))
    print("\n")