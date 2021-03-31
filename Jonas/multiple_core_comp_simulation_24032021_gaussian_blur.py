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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import hyperspy.api as hs
from coreshellp import *
from specerr import *
from specMapDiff import *
import numpy as np
from loadassign import checkLoadFit
from scipy.ndimage import gaussian_filter
sAgPure = hs.load("./Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCuPure = hs.load("./Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")
sCBack=hs.load("./Spectra/Carbonbackground.msa", signal_type="EDS_TEM")#10 nm film
corrat=np.linspace(0,1,11);
ratios=np.linspace(0,1,11);
errt=[];
cerrt=[];
serrt=[];
kfacs = [1,0.72980399]
coreCQ=[];
shellCQ=[];
dens = 20**-1
thickness=0
dim=2
for k in corrat:
    cores=np.zeros((1,11,2048));
    shells=np.zeros((1,11,2048))
    
    for i in range(len(ratios)):
        cores[0][i]=k*sCuPure.data+(1-k)*sAgPure.data
        shells[0][i]=(1-ratios[i])*sAgPure.data+ratios[i]*sCuPure.data
    x=[];
    a=[];
    
    for i in range(len(ratios)):
        x.append(CoreShellP(50,20.0,15.0,dens,dens,1))
        a.append(CoreShellSpec(x[i],cores[0][i],shells[0][i],True))
        a[i].add_background(sCBack, thickness)
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
        a=a.map(gaussian_filter,sigma=2.5)
        cal = hs.load("./Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")#Kalibreringsdata
# For nicer plots, HyperSpy needs some meta data:
    for a in plist:
        a=setCalibration(a, cal)
        cut_spectrum_bottom(a,1000.0)
#Make image
    #imList=[y.get_lines_intensity() for y in plist]
    #for im in imList:
        #redBlueMap(im)
#%%Köra NMF + specAbsErr2D på alla bilder
    #Dessa kommer innehålla felen i total, core och shell:
    err=[]
    coreerr=[]
    shellerr=[]
    #Dessa behövs för kvantifiering senare:
    cFac=[];
    sFac=[];
    for i in range(len(plist)):
        plist[i].decomposition(True,algorithm='NMF',output_dimension =dim) # The "True" variable tells the function to normalize poissonian noise.
        factors = plist[i].get_decomposition_factors() 
        loadings =plist[i].get_decomposition_loadings()
        c,s=checkLoadFit(core[i],shell[i],factors,loadings, dim)
        cFac.append(factors.inav[c])
        sFac.append(factors.inav[s])
        NMFspec1 = cLoadsFacs(loadings, factors)
        NMFparticle1 = NMFspec1.inav[0] + NMFspec1.inav[1]
        err.append(SpecErrAbs2D(NMFparticle1,plist[i]))
        coreerr.append(SpecErrAbs2D(NMFspec1.inav[c],core[i]))
        shellerr.append(SpecErrAbs2D(NMFspec1.inav[s],shell[i]))
    #Sparar felen till totaldata senare.
    cerrt.append(coreerr);
    serrt.append(shellerr);
    errt.append(err);
#%%Kvantifiering
    csF=[hs.stack(cFac),hs.stack(sFac)]
    intensities=[]
    for i in range(2):
        csF[i].set_signal_type("EDS_TEM")
        csF[i].get_calibration_from(cal)
        csF[i].add_elements(['Cu','Ag'])
        csF[i].add_lines()
        bg = csF[i].estimate_background_windows(line_width=[5.0, 7.0])
        intensities.append(csF[i].get_lines_intensity(background_windows=bg))
    
    quant=[]    
    for i in range(len(csF)):
        quant.append(csF[i].quantification(intensities[i], method='CL', factors=kfacs,composition_units='weight'))
    coreCQ.append(quant[0][1].data)
    shellCQ.append(quant[1][1].data)
#%%Plotta
carray=np.array(cerrt)
sarray=np.array(serrt)
earray=np.array(errt)
ccQarr=np.array(coreCQ);
scQarr=np.array(shellCQ);

for i in range(len(ratios)):
    plt.figure(i)
    plt.plot(ratios, carray[i])
    nr=str(10*i)
    plt.xlabel('Fraction Cu in shell')
    plt.ylabel('Relative error in core')
    plt.title(nr+'% Cu in core')
    plt.ylim(0,2)
#Kontur för kärnfelet
K, I = np.meshgrid(100*corrat, 100*ratios)
fig1,ax1=plt.subplots(1,1)
cp1 = ax1.contourf(K, I, carray)
fig1.colorbar(cp1) # Add a colorbar to a plot
ax1.set_title('Relative error (core)')
ax1.set_xlabel('% Cu (shell)')
ax1.set_ylabel('% Cu (core)')#Varje lista i carray hamnar på ett värde på y
plt.show()
#Kontur för skalfelet
fig2,ax2=plt.subplots(1,1)
cp2=ax2.contour(K,I,sarray)
fig2.colorbar(cp2)
ax2.set_title('Relative error (shell)')
ax2.set_xlabel('% Cu (shell)')
ax2.set_ylabel('% Cu (core)')#Varje lista i sarray hamnar på ett värde på y
plt.show()
#Kontur för bedömd kopparhalt (kärna)
fig3,ax3=plt.subplots(1,1)
cp3=ax3.contour(K,I,ccQarr)
fig3.colorbar(cp2)
ax3.set_title('Estimated copper content(core)')
ax3.set_xlabel('% Cu (shell)')
ax3.set_ylabel('% Cu (core)')#Varje lista i sarray hamnar på ett värde på y
plt.show()
#Kontur för bedömd kopparhalt (skal)
fig4,ax4=plt.subplots(1,1)
cp4=ax4.contour(K,I,scQarr)
fig4.colorbar(cp2)
ax4.set_title('Estimated copper content(shell)')
ax4.set_xlabel('% Cu (shell)')
ax4.set_ylabel('% Cu (core)')#Varje lista i sarray hamnar på ett värde på y
plt.show()