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

def CompactRot(factors, loadings, bound_interval):
    angle=RotAngle(loadings.data, bound_interval)
    new_load0=np.cos(angle)*loadings.data[0]+np.sin(angle)*loadings.data[1]
    new_load1=np.cos(angle)*loadings.data[1]-np.sin(angle)*loadings.data[0]
    new_fac0=np.cos(angle)*factors.data[0]+np.sin(angle)*factors.data[1]
    new_fac1=np.cos(angle)*factors.data[1]-np.sin(angle)*factors.data[0]
    loadings.data=np.dstack(new_load0,new_load1)
    factors.data=np.vstack(new_fac0,new_fac1)
    return factors,loadings
    