# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 14:38:15 2021

@author: Jonas
"""

import hyperspy.api as hs
from specMapDiff import *
from specerr import *
import numpy as np

def setCalibration(ucMap,calSpec):
    ucMap.set_signal_type("EDS_TEM")
    ucMap.axes_manager[0].name = 'y'
    ucMap.axes_manager[1].name = 'x'
    ucMap.axes_manager['x'].units = 'nm'
    ucMap.axes_manager['y'].units = 'nm'
    ucMap.axes_manager[-1].name = 'E'
    ucMap.get_calibration_from(calSpec)
    ucMap.add_elements(['Ag','Cu', 'C'])
    ucMap.add_lines(['Ag_La','Cu_Ka','C_Ka'])
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
    combined=combined.transpose(signal_axes=[0],navigation_axes=[1,2,3])
    s = loads.metadata.General.title.split("of")
    combined.metadata.General.title = "Reconstructed: " + s[1]
    return combined

def checkLoadOnly(core,shell,statload,components=2,method='abs'):
   statcoretest = None
   statshelltest = None
   statcore = 1000
   statshell = 1000
   if method=='abs':
       for i in range(components):
           statcoretest=SpecErrAbs2D(statload.inav[i].data, core,signal=False)
           statshelltest=SpecErrAbs2D(statload.inav[i].data, shell,signal=False)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
                statshell = statshelltest
                index_s = i
           if index_c == index_s:
               print('Warning! It guessed the same component for core and shell.')
   elif method=='neuc':
       print("neuc hasn't been implemented for this function yet!")
       for i in range(components):
           statcoretest=SpecErrNEuc(statload.inav[i], core)
           statshelltest=SpecErrNEuc(statload.inav[i], shell)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
               statshell = statshelltest
               index_s = i
           if index_c == index_s:
              print('Warning! It guessed the same component for core and shell.')
   else:
       print(method+" is not a valid comparison method. Enter 'abs' or 'neuc'")
   return index_c,index_s #index för den som bäst passar core respektive shell

def checkLoadFit(core,shell,statfac,statload,components=2,method='abs'):
    
   '''Takes the core and shell of a core-shell particle (Hyperspy-signals) as well as \n
   the factors and loadings of a statistical method and returns the indices\n
   of the factor+loading combo that best correspond to core and shell respectively.
    '''
   Statspec = cLoadsFacs(statload, statfac)
   statcoretest = None
   statshelltest = None
   statcore = np.inf
   statshell = np.inf
   if method=='abs':
       for i in range(components):
           statcoretest=SpecErrAbs2D(Statspec.inav[i], core)
           statshelltest=SpecErrAbs2D(Statspec.inav[i], shell)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
                statshell = statshelltest
                index_s = i
           if index_c == index_s:
               print('Warning! It guessed the same component for core and shell.')
   elif method=='neuc':
        for i in range(components):
           statcoretest=SpecErrNEuc(Statspec.inav[i], core)
           statshelltest=SpecErrNEuc(Statspec.inav[i], shell)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
               statshell = statshelltest
               index_s = i
        if index_c == index_s:
              print('Warning! It guessed the same component for core and shell.')
   else:
       print(method+" is not a valid comparison method. Enter 'abs' or 'neuc'")
   return index_c,index_s #index för den som bäst passar core respektive shell
