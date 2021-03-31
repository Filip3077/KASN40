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

def orBlueMapCuAg2(factors,loadings, title1,title2):
    hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
    plt.title(title1)
    plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
    plt.axvline(8040, c='k', ls=':', lw=0.5)
    plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
    plt.axvline(930, c='k', ls=':', lw=0.5)
    plt.axvline(2984, c='k', ls=':', lw=0.5)
    plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    hs.plot.plot_images(loadings,suptitle=title2, cmap='mpl_colors',
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

        