# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:39:07 2021

@author: Filip
"""
import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt

def specMapDiff(map1,refMap):
    #Om map1 och map2 är hyperspy objekt går det helt enkelt att ta differensen direkt samt att ta absolutvärdet av denna. Om dimentionerna stämmer dvs. 
    diff = abs(map1-refMap)
    return diff


def rel(EF_map,ref):
    '''
    Tar den relativa 
    '''
    refMap = EF_map.data/ref.data
    where_are_NaNs = np.isnan(refMap)
    refMap[where_are_NaNs] = 0
    refMap = hs.signals.BaseSignal(refMap).T
    return refMap




        
