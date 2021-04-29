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
from coreshellp import *
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
cal = hs.load("./Spectra/20 nm cube Fe SSD.msa",signal_type="EDS_TEM")#Kalibreringsdata
tquant=[]
true_tquant=[]
terr=[]
tcerr=[]
tserr=[]
rep=10;
#sFe=hs.load("./Spectra/20 nm cube Fe SSD.msa", signal_type="EDS_TEM")
#sAu=hs.load("./Spectra/20 nm cube Au.msa", signal_type="EDS_TEm")
#cal = hs.load("./Spectra/20nm cube Cu20Ag80.msa",signal_type="EDS_TEM")
#sAgPure=setCalibration(sAgPure,cal)
#sCuPure=setCalibration(sCuPure,cal)
#sCBack=setCalibration(sCBack,cal)
k=0.01*float(input("Input core copper fraction (%):"))
dens = 20**-1
thickness=2
dim=int(input("Input decomposition dimension :"))
gb=bool(int(input("Gaussian blur? (1/0) :")))
tf=bool(int(input("Use radial profile to transfer elements? (1/0):")))
L=len(sAgPure.inav)
ratios=np.linspace(0,1,11);
cores=np.zeros((11,L));
shells=np.zeros((11,L))



for i in range(len(ratios)):#For some reason range(ratios) does not work
    cores[i]=k*sCuPure.data+(1-k)*sAgPure.data
    shells[i]=(1-ratios[i])*sAgPure.data+ratios[i]*sCuPure.data
core_sig=[]
shell_sig=[]
for i in range(len(ratios)):
    cspec=hs.signals.Signal1D(cores[i])
    sspec=hs.signals.Signal1D(shells[i])
    cspec.set_signal_type("EDS_TEM")
    cspec.get_calibration_from(cal)
    cspec.add_elements(['Cu','Ag'])
    cspec.add_lines()
    sspec.set_signal_type("EDS_TEM")
    sspec.get_calibration_from(cal)
    sspec.add_elements(['Cu','Ag'])
    sspec.add_lines()
    core_sig.append(cspec)
    shell_sig.append(sspec)
#%%
for u in range(rep):
    x=[];
    a=[];
    for i in range(len(ratios)):
        x.append(CoreShellP(50,20.0,15.0,dens,dens,1))
        a.append(CoreShellSpec(x[i],cores[i],shells[i],False))
        a[i].add_background(sCBack,thickness)
    #a[i].add_core_component(sFe,0.5)
    #a[i].add_shell_component(sAu,0.1)
# CoreShellP generates two 3D matrices of a sphere. One consisting of the core and one as the shell. 
 # The density here can be seen as having the unit nm^-3 to make the values in the sphere matrix unitless.
# 50x50 pixels, 20nm outer radius, 15nm core radius, densities, 1x1 nm pixel size.

# CoreShellSpec fills the core and shell matrices from above with the simulated spectra and are then turned into 
# HyperSpy objects for the latter comparrisons. Combining these gives us a HyperSpy object of a whole particle (p). 
    core = [y.core for y in a]
    shell = [y.shell for y in a]
    bcore=[y.base.core for y in a]
    bshell=[y.base.shell for y in a]
#core=list(map(lambda x: hs.signals.Signal1D(x),core))
#shell=list(map(lambda x: hs.signals.Signal1D(x),shell))
    parts=[y.getmatr() for y in a]
    plist=list(map(lambda x: hs.signals.Signal1D(x),[x.data for x in parts]));
    clist=list(map(lambda x: hs.signals.Signal1D(x),core))
    slist=list(map(lambda x: hs.signals.Signal1D(x), shell))
#plist=parts

    for a in plist:
        a.add_poissonian_noise()# Adds poissonian noise to the existing spectra.
        if gb==True:
            a=a.map(gaussian_filter,sigma=2)
        

# For nicer plots, HyperSpy needs some meta data:
    for a in plist:
        a=setCalibration(a, cal)
        cut_spectrum_bottom(a,1500.0)
        cut_spectrum_top(a,8300.0)
        a=a.rebin([2,2,10])
#Make image
    imList=[y.get_lines_intensity() for y in plist]
#for im in imList:
    #redBlueMap(im)
#%%Köra NMF + specAbsErr2D på alla bilder
    err=[]
    coreerr=[]
    shellerr=[]
    cFac=[];
    sFac=[];
    NMFparts=[];
    cchoice=[];
    schoice=[];
    kfacs = [1,0.72980399]
    for i in range(len(plist)):
        plist[i].decomposition(True,algorithm='NMF',output_dimension =dim) # The "True" variable tells the function to normalize poissonian noise.
        factors = plist[i].get_decomposition_factors() 
        loadings =plist[i].get_decomposition_loadings()
    #c,s=0,1;
        if tf==True:
            c,s=checkLoadFit(clist[i],slist[i],factors,loadings, dim,'abs')
            factors,loadings=transfer_elements(factors,loadings,c,s,50)
        c,s=checkLoadFit(clist[i],slist[i],factors,loadings, dim,'abs')
        cchoice.append(c);
        schoice.append(s);
        NMFspec1 = cLoadsFacs(loadings, factors)
        cFac.append(factors.inav[c]);
        sFac.append(factors.inav[s]);
        NMFparticle1 = NMFspec1.inav[c] + NMFspec1.inav[s]
        #orBlueMapCuAg(factors,loadings,'NMF')
        err.append(SpecErrAbs2D(NMFparticle1,plist[i]));
        coreerr.append(SpecErrAbs2D(NMFspec1.inav[c],clist[i]));
        shellerr.append(SpecErrAbs2D(NMFspec1.inav[s],slist[i]));
        NMFparts.append(NMFparticle1);
    
#%%Kvantifiera
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
    tquant.append(quant)
    terr.append(err)
    tcerr.append(coreerr)
    tserr.append(shellerr)
#%%Kvantifiera verkliga partiklar
true_csF=[hs.stack(core_sig), hs.stack(shell_sig)];
intensities=[];
for i in range(2):
    true_csF[i].set_signal_type("EDS_TEM")
    true_csF[i].get_calibration_from(cal)
    true_csF[i].add_elements(['Cu','Ag'])
    true_csF[i].add_lines()
    true_bg=true_csF[i].estimate_background_windows(line_width=[5.0,7.0])
    intensities.append(true_csF[i].get_lines_intensity(background_windows=bg))
true_quant=[]
for i in range(len(true_csF)):
    true_quant.append(true_csF[i].quantification(intensities[i],method='CL',factors=kfacs,composition_units='weight'))
true_tquant.append(true_quant)
#%%Statistik för Cu-innehåll

core_cu=np.zeros((rep,11))
shell_cu=np.zeros((rep,11))

for i in range(len(tquant)):
    for q in range(len(tquant[0][0][1].data)):
        core_cu[i][q]+=tquant[i][0][1].data[q]
        shell_cu[i][q]+=tquant[i][1][1].data[q]
mean_core_cu=np.sum((1/rep)*core_cu,axis=0);
mean_shell_cu=np.sum((1/rep)*shell_cu,axis=0)
core_std=np.std(core_cu,axis=0);
shell_std=np.std(shell_cu,axis=0);

#%%Statistik för absolutfel
core_err=np.zeros((rep,11))
shell_err=np.zeros((rep,11))
tot_err=np.zeros((rep,11))

for i in range(len(terr)):
    for q in range(len(terr[0])):
        core_err[i][q]+=tcerr[i][q]
        shell_err[i][q]+=tserr[i][q]
        tot_err[i][q]+=terr[i][q]
mean_core_err=np.sum((1/rep)*core_err,axis=0);
mean_shell_err=np.sum((1/rep)*shell_err,axis=0);
mean_err=np.sum((1/rep)*tot_err,axis=0);
core_err_std=np.std(core_err,axis=0);
shell_err_std=np.std(shell_err,axis=0);
tot_err_std=np.std(tot_err,axis=0);

#%%Plotta
np.array(coreerr)
plt.figure(1001)
plt.errorbar(100*ratios, mean_core_err,1.96*core_err_std*1/np.sqrt(rep),fmt='.k')
plt.plot(100*ratios,mean_core_err)
plt.ylim(0,2)
cont=str(100*k)
plt.title('Relative error of the core with '+cont+'% Cu in core')
plt.xlabel('Percentage Cu in shell')
plt.ylabel('Relative error (core)')
plt.figure(1002)
plt.plot(ratios,cchoice)
plt.figure(1003)
plt.errorbar(100*ratios, mean_core_cu,1.96*core_std*1/np.sqrt(rep),fmt='.k');
plt.plot(100*ratios,mean_core_cu);
plt.plot(100*ratios,np.ones(11)*k*100)
plt.xlabel('Percentage Cu in shell')
plt.ylabel('CL estimate of core Cu content(%)')
plt.title('CL estimate of core Cu content with a true value of '+cont+'%')
plt.ylim([0, 100])
plt.figure(1004)
plt.errorbar(100*ratios,mean_shell_cu,1.96*shell_std*1/np.sqrt(rep),fmt='.k')
plt.plot(100*ratios, mean_shell_cu)
plt.plot(100*ratios,ratios*100)
plt.xlabel('Percentage Cu in shell')
plt.ylabel('CL estimate of shell Cu content(%)')
plt.title('CL estimate of shell Cu content with a true core Cu content of ' +cont+'%')