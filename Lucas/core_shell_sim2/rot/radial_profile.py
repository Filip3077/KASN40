# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 16:49:39 2021

@author: Jonas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

def CorrectShell(b, *args):
    temp_shell = args[1]+b*args[0]
    temp_profile = radial_profile(temp_shell, [24.5, 24.5])
    
    return np.sum((temp_profile/np.max(temp_profile)-args[2])**2)

def transfer_elements(factors,loadings,c,s,size):
    dest_tot = np.sqrt(np.sum((loadings.inav[c].data+loadings.inav[s].data >0.01), axis=(0,1))/np.pi)
    dest_core = np.sqrt(np.sum((loadings.inav[c].data >0.01), axis=(0,1))/np.pi)
    dest_shell = dest_tot-dest_core

    center_int_frac = dest_shell/np.sqrt(dest_tot**2-dest_core**2)

    shell_profile = radial_profile(loadings.inav[1].data, [(size-1)/2, (size-1)/2])
    shell_profile_ideal = np.arange(0,shell_profile.shape[0],1)**2
    shell_profile_ideal = np.nan_to_num(np.sqrt(dest_tot**2-shell_profile_ideal)) - np.nan_to_num(np.sqrt(dest_core**2-shell_profile_ideal))
    shell_profile_ideal= shell_profile_ideal / np.max(shell_profile_ideal)
    loadings_tup=[loadings.inav[c].data,loadings.inav[s].data]
    res1 = minimize_scalar(CorrectShell, args=tuple(loadings_tup)+(shell_profile_ideal,), bounds=(-0.8, 0.8), method='bounded')
    b = res1.x
    loadings_corr=loadings
    factors_corr=factors
    loadings_corr.inav[s].data=loadings.inav[s].data+b*loadings.inav[c].data
    factors_corr.inav[c].data=factors.inav[c].data-b*factors.inav[s].data
    return factors_corr,loadings_corr
    