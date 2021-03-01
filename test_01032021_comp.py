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

sAgPure = hs.load("./Spectra/20nm cube Cu0Ag100.msa", signal_type="EDS_TEM")
sCuPure = hs.load("./Spectra/20nm cube Cu100Ag0.msa", signal_type="EDS_TEM")

ratios = np.linspace(0, 1, 11);
cores = np.zeros((11, 11, 2048));
shells = np.zeros((11, 11, 2048))

for i in range(len(ratios)):  # For some reason range(ratios) does not work
    for j in range(len(ratios)):
        cores[i][j][0:2048] = (1-ratios[i])*sCuPure.data+ratios[i]*sAgPure.data
        shells[i][j][0:2048] = (1-ratios[j])*sAgPure.data+ratios[j]*sCuPure.data
x = [];
a = [];
dens = 20**-1
for i in range(len(ratios)):
    x.append([])
    a.append([])
    for j in range(len(ratios)):
        x[i].append(CoreShellP(50, 20.0, 15.0, dens, dens, 1))
        a[i].append(CoreShellSpec(x[i][j], cores[0][i][j], shells[0][i][j], False))
# CoreShellP generates two 3D matrices of a sphere. One consisting of the core and one as the shell.
 # The density here can be seen as having the unit nm^-3 to make the values in the sphere matrix unitless.
# 50x50 pixels, 20nm outer radius, 15nm core radius, densities, 1x1 nm pixel size.

# CoreShellSpec fills the core and shell matrices from above with the simulated spectra and are then turned into
# HyperSpy objects for the latter comparrisons. Combining these gives us a HyperSpy object of a whole particle (p).
core=[] [y.core for y in a]
shell=[] [y.shell for y in a]
parts=[]
for i in range(len(ratios)):
    core.append([])
    shell.append([])
    parts.append([])
    for j in range(len(ratios)):
        core[i].append(a[i][j].core)
        shell[i].append(a[i][j].shell)
        parts.append(a[i][j].getmatr())
# core=list(map(lambda x: hs.signals.Signal1D(x),core))
# shell=list(map(lambda x: hs.signals.Signal1D(x),shell))
parts = [y.getmatr() for y in a]
plist= list(map(lambda x: hs.signals.Signal1D(x), [x.data for x in parts]));
clist=list(map(lambda x: hs.signals.Signal1D(x), core))
slist=list(map(lambda x: hs.signals.Signal1D(x), shell))
# plist=parts
for a in plist:
    a.add_poissonian_noise()  # Adds poissonian noise to the existing spectra.
cal=hs.load("./Spectra/20nm cube Cu20Ag80.msa",
            signal_type="EDS_TEM")  # Kalibreringsdata
# For nicer plots, HyperSpy needs some meta data:
for a in plist:
    a=setCalibration(a, cal)
# Make image
imList=[y.get_lines_intensity() for y in plist]
for im in imList:
    redBlueMap(im)
# %%Köra NMF + specAbsErr2D på alla bilder
err=[]
coreerr=[]
shellerr=[]
dim=2
NMFparts=[];
for i in range(len(plist)):
    # The "True" variable tells the function to normalize poissonian noise.
    plist[i].decomposition(True, algorithm='NMF', output_dimension=dim)
    factors=plist[i].get_decomposition_factors()
    loadings=plist[i].get_decomposition_loadings()
    NMFspec1=cLoadsFacs(loadings, factors)
    NMFparticle1=NMFspec1.inav[0] + NMFspec1.inav[1]
    orBlueMapCuAg(factors, loadings, 'NMF')
    err.append(SpecErrAbs2D(NMFparticle1, plist[i]))
    coreerr.append(SpecErrAbs2D(NMFspec1.inav[0], core[i]))
    shellerr.append(SpecErrAbs2D(NMFspec1.inav[1], shell[i]))
    NMFparts.append(NMFparticle1);
