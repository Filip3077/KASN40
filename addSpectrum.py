# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 18:58:23 2021

@author: Filip
"""
import hyperspy.api as hs
import numpy as np

def addSpectrum(a,spec,specthickness):
    L = len(spec.data)
    mat = np.zeros((len(a),len(a[0]),L))

    for i in range(0,len(a)):
        for j in range(0,len(a[0])):
            mat[i,j,0:L]=(a[i,j]/specthickness)*spec.data
            
    specMat = hs.signals.Signal1D(mat)
    specMat.set_signal_type("EDS_TEM")
    specMat.get_calibration_from(spec)
    return specMat