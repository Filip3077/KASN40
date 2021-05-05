# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 15:31:58 2021

@author: Lucas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from coreshellfunctions import*
#%%

sAg = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCu = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")
sC=hs.load("../Spectra/Carbonbackground.msa", signal_type="EDS_TEM")
cal = sAg #0.45*sAg.data+0.45*sCu.data+0.1*sC.data

rep = 10
size = 50
intensity = 1/20
spectrum_cutoff = 500
a = CoreShellP(size, 20, 15, intensity, intensity, 1)
newbin = [1, 1, 8]
ps = []; cores = []; shells =[]
for i in range(10,49,5):
    specCore = sCu.data*(1-i/100) + sAg.data*i/100
    specShell = sAg.data*(1-i/100) + sCu.data*i/100
    b = CoreShellSpec(a, specCore,specShell)
    b.add_background(sC,1)
    p_temp = b.core + b.shell
    core_temp = b.core
    shell_temp = b.shell
    p_temp.add_poissonian_noise()
    p_temp = setCalibration(p_temp, cal); core_temp = setCalibration(core_temp, cal); shell_temp = setCalibration(shell_temp, cal)
    cut_spectrum_bottom(p_temp, spectrum_cutoff); cut_spectrum_bottom(core_temp, spectrum_cutoff); cut_spectrum_bottom(shell_temp, spectrum_cutoff)
    p_temp.change_dtype('float64'); core_temp.change_dtype('float64'); shell_temp.change_dtype('float64')
    p_temp = p_temp.rebin(scale = newbin); core_temp = core_temp.rebin(scale = newbin); shell_temp = shell_temp.rebin(scale = newbin)
    ps.append(p_temp); cores.append(core_temp); shells.append(shell_temp)
    
#%%

nfac = 2 #int(input("Chose nr of components :"))

loadings = []; factors = []; factors_corr = []; loadings_corr = []
NMF_ps = []; NMF_corr_ps =[]
wtp = [[],[]]; wtp_corr = [[],[]]
f = 0

for i in range(len(ps)):
    print('PARTICLE: '+str(i))
    for r in range(rep):
        print('REPETITION:'+str(f+r))
        ps[i].decomposition(True, algorithm='NMF', output_dimension=nfac)
        loadings.append(ps[i].get_decomposition_loadings())
        factors.append(ps[i].get_decomposition_factors())
        factors_corr_temp, loadings_corr_temp = transfer_elements(factors[f+r], loadings[f+r], 0,1)
        factors_corr.append(factors_corr_temp); loadings_corr.append(loadings_corr_temp)
        
        NMF_p = cLoadsFacs(loadings[f+r], factors[f+r]).split()
        NMF_corr_p = cLoadsFacs(loadings_corr[f+r], factors_corr[f+r]).split()
        for j in range(2):
            NMF_p[j] = setCalibration(NMF_p[j], ps[i].inav[24,24] )
            NMF_corr_p[j] = setCalibration(NMF_corr_p[j], ps[i].inav[24,24])
        NMF_ps.append(NMF_p); NMF_corr_ps.append(NMF_corr_p)
        
        wtp[0].append(quantify(NMF_p[0].inav[:,:].sum())); wtp[1].append(quantify(NMF_p[1].inav[:,:].sum()))
        wtp_corr[0].append(quantify(NMF_corr_p[0].inav[:,:].sum())); wtp_corr[1].append(quantify(NMF_corr_p[1].inav[:,:].sum()))
        
    meanAgCore.append(np.sum([t[0] for t in wtp[0][f:f+rep]])/rep)
    meanCuCore.append(np.sum([t[1] for t in wtp[0][f:f+rep]])/rep)
    meanAgShell.append(np.sum([t[0] for t in wtp[1][f:f+rep]])/rep)
    meanCuShell.append(np.sum([t[1] for t in wtp[0][f:f+rep]])/rep)
    steAgCore.append(np.std([t[0] for t in wtp[0][f:f+rep]], ddof=1)/np.sqrt(np.size([t[0] for t in wtp[0][f:f+rep]]))    
    f+=rep

# wtp_means = []
# for i in range(len(ps)):
#     meanAg = np.sum(wtp[0][i:i+9])/rep
#     meanCu = np.sum(wtp[1][i:i+9])/rep
#     wtp_means.append([meanAg, meanCu]) 
            
    
#%%


    
