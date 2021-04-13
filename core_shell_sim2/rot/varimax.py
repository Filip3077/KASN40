# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:39:32 2020

@author: Martin
"""


import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

def varimax(Phi , gamma = 1.0 , q = 100  , tol = 1e-10):
    ''' The original Varimax function given by Martin\n
    Takes the unfolded (2+1D -> 1+1D) loadings matrix as main variable. '''
    import numpy as np
    from numpy import linalg
    p, k = Phi.shape
    R = np.eye(k)
    d=1e-6
    
    for i in range(q) :
        d_old = d
        Lambda = np.dot( Phi , R)
        temp1 = np.dot( Lambda , np.diag( np.diag( np.dot( Lambda.T , Lambda ) ) ) )
        temp2 = np.dot( Phi.T , np.asarray( Lambda )**3 - ( gamma/ p ) * temp1 )
        u,s,vh = linalg.svd(temp2)
        R = np.dot ( u , vh )
        d = np.sum( s )
        
        # print(d/d_old)
        
        if d/d_old < (1+tol) :
            break
    return R

def FullVarimax(factors,loadings,components=-1):
    ''' Takes factors and loadings and returns (positive) rotated factors and loadings\n
    Optional variable components should be used if the number of components have not been choosen before.'''
    import numpy as np
    from numpy import linalg
    
    if components == -1:
        nfac = len(factors)
    else:
        nfac = components
        
    factors_selected = factors.inav[0:nfac]
    loadings_selected = loadings.inav[0:nfac]

    #Unfold to turn loadings 2+1D matrix into 1+1D matrix for varimaxfunction.
    loadings_unfold = loadings_selected.deepcopy()
    loadings_unfold.unfold()
    factors_temp = factors_selected.deepcopy()
    factors_selected = factors_selected.data
    loadings_selected = loadings_unfold.data

    R =varimax(loadings_selected.T)

    loadings_selected_rot = np.matmul(loadings_selected.T, R).T
    factors_selected_rot = np.matmul(linalg.inv(R), factors_selected)
    
    #Flipping negative factors and loadings. Only guaranteed for fully negative (non-mixed negative and positive peaks)
    for i in range(nfac):
        if factors_selected_rot[i,:].sum() < 0:
            factors_selected_rot[i,:] = -factors_selected_rot[i,:]
            if loadings_selected_rot[i,:].sum() < 0:
                loadings_selected_rot[i,:] = -loadings_selected_rot[i,:]
    
    factors_temp.data = factors_selected_rot # this way factors is kept as hs objekt.
    #Fold the loadings matrix to original format.
    loadings_unfold.data = loadings_selected_rot
    loadings_unfold.fold()
    
    return factors_temp, loadings_unfold