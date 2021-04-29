# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 09:15:41 2021

@author: Filip
"""

from core_shell_sim.sim.coreshellp_3D import *
import hyperspy.api as hs
 
s = hs.load("./Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]

f1 =sCu*0.9 + sAg*0.1
f2 = sCu*0.1 + sAg*0.9

mat = CoreShellP_3D(50,10,5,1,1,1)
cs = CoreShellSpec_3D(mat,f1,f2)
s = cs.getmatr()
s.add_poissonian_noise(keep_dtype=True)


s.decomposition()
s.plot_explained_variance_ratio()

s.decomposition(algorithm = 'NMF', output_dimension = 2)
s_loads = s.get_decomposition_loadings()
s_facs = s.get_decomposition_factors()


hs.plot.plot_spectra(s_facs,style='cascade',padding=-1)

l1 = s_loads.inav[0].T
l2 = s_loads.inav[1].T

vis1 = [l1.inav[25,:,:],l1.inav[20,:,:],l1.inav[18,:,:]]
vis2 = [l2.inav[25,:,:],l2.inav[20,:,:],l2.inav[18,:,:]]

hs.plot.plot_images(vis1, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})

hs.plot.plot_images(vis2, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})