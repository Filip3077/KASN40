# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:39:07 2021

@author: Filip
"""
import hyperspy.api as hs
import numpy as np

def specMapDiff(map1,map2):
    '''Generates a "difference map" from hyperspy signal objects map1 and map2'''
    #Om map1 och map2 är hyperspy objekt går det helt enkelt att ta differensen direkt samt att ta absolutvärdet av denna. Om dimentionerna stämmer dvs. 
    diff = abs(map1-map2)
    return diff
    

def rel(EF_map,ref):
    '''EF_map is assumed to be a difference map of intensities as a hyperspy image\n
    i.e. simulation-ref or similar. ref is the reference image as a hyperspy image'''
    #100:100:1
    size  = len(EF_map.inav[0])
    for x in range(size):
        for y in range(size):
            if (ref.inav[x,y].data == 0):
                EF_map.inav[x,y] = 0
            else:
                EF_map.inav[x,y] = EF_map.inav[x,y].data/ref.inav[x,y].data
                
            
            
    return EF_map

def setCalibration(ucMap,refMap):
    '''Executes a a standard set of hyperspy signal commands for our project'''
    ucMap.set_signal_type("EDS_TEM")
    ucMap.axes_manager[0].name = 'y'
    ucMap.axes_manager[1].name = 'x'
    ucMap.axes_manager['x'].units = 'nm'
    ucMap.axes_manager['y'].units = 'nm'
    ucMap.axes_manager[-1].name = 'E'
    ucMap.get_calibration_from(refMap)
    ucMap.add_elements(['Ag','Cu'])
    ucMap.add_lines(['Ag_La','Cu_Ka'])
    return ucMap

def cLoadsFacs(loads,facs):
    #Antar att både loads och facs kommer från samma "ursprung" och har samma ordning och dimentioner
    #För att få ett korrekta dimentioner på  hyperspy objektet böhöver loads transponeras från [| x y]  till [x y |] har att göra med hur energiaxeln behandlas 
    
    dim = len(loads)
    size = len(loads.isig)
    esize = len(facs.isig)
    combinedMat = np.empty((dim,size,size,esize))
    
    for i in range(dim):
        #För att få ett korrekta dimentioner på  hyperspy objektet behöver loads transponeras från [| x y]  till [x y |] har att göra med hur energiaxeln behandlas 
        combinedMat[i] = (loads.inav[i].T*facs.inav[i]).data
    
    
    combined = hs.signals.BaseSignal(combinedMat)
    combined=combined.transpose(signal_axes=[0],navigation_axes=[3, 2, 1])
    return combined

        