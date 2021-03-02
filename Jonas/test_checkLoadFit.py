# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:17:56 2021

@author: Jonas
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec
from specerr import *
from specMapDiff import *
import numpy as np
from coreshellFunctions import checkLoadFit
sAgPure = hs.load("./Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("./Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")
cal = hs.load("./Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")#
dens=20**-1;
x=CoreShellP(50,20.0,15.0,dens,dens,1);
a=CoreShellSpec(x,sCuPure,sAgPure,True);
p=a.getmatr()
p.add_poissonian_noise();
p=setCalibration(p,cal)
#%%NMF-decomp
p.decomposition(True,algorithm='NMF',output_dimension =2) # The "True" variable tells the function to normalize poissonian noise.
factors = p.get_decomposition_factors() 
loadings =p.get_decomposition_loadings()
    #c,s=0,1;
c,s=checkLoadFit(a.core,a.shell,factors,loadings, 2)
NMFspec1 = cLoadsFacs(loadings, factors)
NMFparticle1 = NMFspec1.inav[0] + NMFspec1.inav[1]
orBlueMapCuAg(factors,loadings,'NMF')
err=SpecErrAbs2D(NMFparticle1,p)
coreerr=SpecErrAbs2D(NMFspec1.inav[c],a.core)
shellerr=SpecErrAbs2D(NMFspec1.inav[s],a.shell)

