# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 10:33:59 2021

@author: Jonas
"""

import numpy as np
import hyperspy.api as hs
from scipy.optimize import minimize_scalar

def RotCompactness(a, *loads):

    return np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1]).sum(axis=(0,1))

def RotAngle(loads,bound_interval):
    res = minimize_scalar(RotCompactness, args=loads, bounds=bound_interval, method='bounded')
    return res.x

def CompactRot(factors, loadings, bound_interval,nfac):
    loads=(loadings.data[0,:,:],loadings.data[1,:,:])
    a1=RotAngle(loads, bound_interval)
    R1 = np.eye(nfac); 
    R1[0,:] = [np.cos(a1), -np.sin(a1)]; R1[1,:] = [np.sin(a1), np.cos(a1)]
    rotload=loadings.deepcopy();
    rotfac=factors.deepcopy();
    rotload.data=np.matmul(loadings.data.T,R1).T
    rotfac.data=np.matmul(factors.data.T,R1).T;
    return rotfac,rotload
    