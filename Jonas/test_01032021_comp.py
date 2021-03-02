# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 14:46:46 2021

@author: Jonas

"""
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 14:46:46 2021

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
k=0.01*float(input("Input core copper fraction (%):"))

ratios=np.linspace(0,1,11);
cores=np.zeros((1,11,2048));
shells=np.zeros((1,11,2048))

for i in range(len(ratios)):#For some reason range(ratios) does not work
    cores[0][i]=k*sCuPure.data+(1-k)*sAgPure.data
    shells[0][i]=(1-ratios[i])*sAgPure.data+ratios[i]*sCuPure.data
x=[];
a=[];
dens = 20**-1
for i in range(len(ratios)):
    x.append(CoreShellP(50,20.0,15.0,dens,dens,1))
    a.append(CoreShellSpec(x[i],cores[0][i],shells[0][i],True))
# CoreShellP generates two 3D matrices of a sphere. One consisting of the core and one as the shell. 
 # The density here can be seen as having the unit nm^-3 to make the values in the sphere matrix unitless.
# 50x50 pixels, 20nm outer radius, 15nm core radius, densities, 1x1 nm pixel size.

# CoreShellSpec fills the core and shell matrices from above with the simulated spectra and are then turned into 
# HyperSpy objects for the latter comparrisons. Combining these gives us a HyperSpy object of a whole particle (p). 
core = [y.core for y in a]
shell = [y.shell for y in a]
#core=list(map(lambda x: hs.signals.Signal1D(x),core))
#shell=list(map(lambda x: hs.signals.Signal1D(x),shell))
parts=[y.getmatr() for y in a]
plist=list(map(lambda x: hs.signals.Signal1D(x),[x.data for x in parts]));
clist=list(map(lambda x: hs.signals.Signal1D(x),core))
slist=list(map(lambda x: hs.signals.Signal1D(x), shell))
#plist=parts
for a in plist:
    a.add_poissonian_noise()# Adds poissonian noise to the existing spectra.
cal = hs.load("./Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")#Kalibreringsdata
# For nicer plots, HyperSpy needs some meta data:
for a in plist:
    a=setCalibration(a, cal)
#Make image
imList=[y.get_lines_intensity() for y in plist]
for im in imList:
    redBlueMap(im)
#%%Köra NMF + specAbsErr2D på alla bilder
err=[]
coreerr=[]
shellerr=[]
dim=2
NMFparts=[];
for i in range(len(plist)):
    plist[i].decomposition(True,algorithm='NMF',output_dimension =dim) # The "True" variable tells the function to normalize poissonian noise.
    factors = plist[i].get_decomposition_factors() 
    loadings =plist[i].get_decomposition_loadings()
    #c,s=0,1;
    c,s=checkLoadFit(core[i],shell[i],factors,loadings, dim)
    NMFspec1 = cLoadsFacs(loadings, factors)
    NMFparticle1 = NMFspec1.inav[0] + NMFspec1.inav[1]
    orBlueMapCuAg(factors,loadings,'NMF')
    err.append(SpecErrAbs2D(NMFparticle1,plist[i]))
    coreerr.append(SpecErrAbs2D(NMFspec1.inav[c],core[i]))
    shellerr.append(SpecErrAbs2D(NMFspec1.inav[s],shell[i]))
    NMFparts.append(NMFparticle1);
#%%Plotta
np.array(coreerr)
plt.plot(ratios, coreerr)