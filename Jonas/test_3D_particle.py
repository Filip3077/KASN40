# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 16:13:37 2021

@author: Jonas
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellp import *
from coreshellp_3D import *
from specerr import *
from specMapDiff import *
import numpy as np
from loadassign import checkLoadFit
from scipy.ndimage import gaussian_filter
from k_factors import silver_k_factor
from radial_profile import transfer_elements

sAgPure = hs.load("./Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("./Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")
sCBack=hs.load("./Spectra/Carbonbackground.msa", signal_type="EDS_TEM")
cal = hs.load("./Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")#hs.load("./Spectra/20 nm cube Fe SSD.msa",signal_type="EDS_TEM")
kfacs = [1,0.72980399]
x=CoreShellP_3D(50,20.0,15.0,1,1,1);
core_spec=0.1*sCuPure.data+0.9*sAgPure.data
shell_spec=0.1*sAgPure.data+0.9*sCuPure.data
a=CoreShellSpec_3D(x,core_spec,shell_spec,False)
b=a.getmatr()
p=hs.signals.Signal1D(b.data)

p.set_signal_type("EDS_TEM")
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager[2].name = 'z'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager['z'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.get_calibration_from(cal)
p.add_elements(['Ag','Cu'])
p.add_lines(['Ag_La','Cu_Ka'])
p=p.map(gaussian_filter,sigma=2.0)

p.decomposition(True, algorithm='NMF',output_dimension=2)
factors=p.get_decomposition_factors()
loadings=p.get_decomposition_loadings()
factors.set_signal_type('EDS_TEM')
factors.get_calibration_from(cal)
factors.add_elements(['Cu','Ag'])
factors.add_lines()
bg = factors.estimate_background_windows(line_width=[5.0, 7.0])
intensities=factors.get_lines_intensity(background_windows=bg)
quant=factors.quantification(intensities, method='CL', factors=kfacs, compositions_units='weight')
#%% Quantification of original particle
core_spec=hs.signals.Signal1D(core_spec)
core_spec.set_signal_type('EDS_TEM')
core_spec.get_calibration_from(cal)
core_spec.add_elements(['Cu','Ag'])
core_spec.add_lines()
bg2=core_spec.estimate_background_windows([5.0,7.0])
intensities2=core_spec.get_lines_intensity(line_width=[5.0, 7.0])
quant2=core_spec.quantification(intensities2, method='CL', factors=kfacs, composition_units='weight')




