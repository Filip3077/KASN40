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
cal = hs.load("./Spectra/20 nm cube Fe SSD.msa",signal_type="EDS_TEM")
kfacs = [1,0.72980399]
x=CoreShellP_3D(50,20.0,10.0,1,1,1);
core_spec=0.9*sCuPure.data+0.1*sAgPure.data
shell_spec=0.9*sAgPure.data+0.1*sCuPure.data
a=CoreShellSpec_3D(x,core_spec,shell_spec,False)
b=a.getmatr()
a=hs.signals.Signal1D(b.data)

a.set_signal_type("EDS_TEM")
a.axes_manager[0].name = 'y'
a.axes_manager[1].name = 'x'
a.axes_manager[2].name = 'z'
a.axes_manager['x'].units = 'nm'
a.axes_manager['y'].units = 'nm'
a.axes_manager['z'].units = 'nm'
a.axes_manager[-1].name = 'E'
a.get_calibration_from(cal)
a.add_elements(['Ag','Cu'])
a.add_lines(['Ag_La','Cu_Ka'])

a.decomposition(True, algorithm='NMF',output_dimension=2)
factors=a.get_decomposition_factors()
loadings=a.get_decomposition_loadings()
factors.set_signal_type('EDS_TEM')
factors.get_calibration_from(cal)
factors.add_elements(['Cu','Ag'])
factors.add_lines()
bg = factors.estimate_background_windows(line_width=[5.0, 7.0])
intensities=factors.get_lines_intensity(background_windows=bg)
quant=factors.quantification(intensities, method='CL', factors=kfacs, compositions_units='weight')



