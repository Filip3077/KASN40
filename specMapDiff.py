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

def redBlueMap(im,supTitle=None, label=None):
    if (supTitle==None and label==None):
        hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
                            colorbar='single', vmin='1th', vmax='99th', scalebar='all',
                            scalebar_color='black', suptitle_fontsize=16,
                            padding={'top':0.8, 'bottom':0.10, 'left':0.05,
                                     'right':0.85, 'wspace':0.20, 'hspace':0.10})
    else:
      hs.plot.plot_images(im, suptitle=supTitle, label=label,tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
                            colorbar='single', vmin='1th', vmax='99th', scalebar='all',
                            scalebar_color='black', suptitle_fontsize=16,
                            padding={'top':0.8, 'bottom':0.10, 'left':0.05,
                                     'right':0.85, 'wspace':0.20, 'hspace':0.10})  
    return None

def orBlueMapCuAg(factors,loadings, title):
    hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
    plt.title(title)
    plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
    plt.axvline(8040, c='k', ls=':', lw=0.5)
    plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
    plt.axvline(930, c='k', ls=':', lw=0.5)
    plt.axvline(2984, c='k', ls=':', lw=0.5)
    plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    hs.plot.plot_images(loadings, cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
    return None

def rel(EF_map,ref):
    '''
    Tar den relativa 
    '''
    refMap = EF_map.data/ref.data
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
    return combined

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
        
