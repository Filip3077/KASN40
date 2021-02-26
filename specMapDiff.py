# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:39:07 2021

@author: Filip
"""
import hyperspy.api as hs
import numpy as np
from specerr import SpecErrAbs2D

def specMapDiff(map1,refMap):
    #Om map1 och map2 är hyperspy objekt går det helt enkelt att ta differensen direkt samt att ta absolutvärdet av denna. Om dimentionerna stämmer dvs. 
    diff = abs(map1-refMap)
    return diff
    

def rel(EF_map,ref):
    '''
    Tar den relativa 
    '''
    refMap = np.divide(EF_map.data, ref.data, out=np.zeros_like(EF_map.data), where=ref.data!=0)
    where_are_NaNs = np.isnan(refMap)
    refMap[where_are_NaNs] = 0
    refMap = hs.signals.BaseSignal(refMap).T
    return refMap


def setCalibration(ucMap,calSpec):
    ucMap.set_signal_type("EDS_TEM")
    ucMap.axes_manager[0].name = 'y'
    ucMap.axes_manager[1].name = 'x'
    ucMap.axes_manager['x'].units = 'nm'
    ucMap.axes_manager['y'].units = 'nm'
    ucMap.axes_manager[-1].name = 'E'
    ucMap.get_calibration_from(calSpec)
    ucMap.add_elements(['Ag','Cu'])
    ucMap.add_lines(['Ag_La','Cu_Ka'])
    return ucMap

def cLoadsFacs(loads,facs):
    #Antar att både loads och facs kommer från samma "ursprung" och har samma ordning och dimentioner
    #För att få ett korrekta dimentioner på  hyperspy objektet böhöver loads transponeras från [| x y]  till [x y |] har att göra med hur energiaxeln behandlas 
    
    dim = len(loads)

    combinedMat = []
    
    for i in range(dim):
        #För att få ett korrekta dimentioner på  hyperspy objektet behöver loads transponeras från [| x y]  till [x y |] har att göra med hur energiaxeln behandlas 
        combinedMat.append(loads.inav[i].T*facs.inav[i])
    return combinedMat

    




    
    
    
    
    
    
    
    
    
    
    