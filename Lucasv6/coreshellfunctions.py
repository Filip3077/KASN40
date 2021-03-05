# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:39:07 2021

@author: Filip
"""
import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from coreshellp import CoreShellP, CoreShellSpec
from specerr import *

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
    combined=combined.transpose(signal_axes=[0],navigation_axes=[3, 2, 1])
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
    loadings_temp = loadings_selected.deepcopy()
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
    
    
    # factors_temp.data = factors_selected_rot # this way factors is kept as a hs objekt.
    #Fold the loadings matrix to original format.
    loadings_unfold.data = loadings_selected_rot
    loadings_unfold.fold()
    
    factors_temp.data = factors_selected_rot
    loadings_temp = loadings_unfold
    return factors_temp, loadings_temp
        
def genfullparticle(size,r1,r2,specCore,specShell,signal=True,dens1=1,dens2=1,l=1):
    ''' Generates a full core@shell particle using coreshellp.py \n
    Only size,r1,r2 and the spectrums is necessary, densities and pixel size \n
    is set to 1 by defult but can be changed.\n
    Whichever of r1 or r2 can be the large or small radius.'''
    a = CoreShellP(size,r1,r2,dens1,dens2,l) 
    return CoreShellSpec(a,specCore,specShell,signal)

def checkLoadFit(self,core,shell,statfac,statload,component=2):
    
   '''Takes the core and shell of a core-shell particle (Hyperspy-signals) as well as \n
   the factors and loadings of a statistical method and returns the indices\n
   of the factor+loading combo that best correspond to core and shell respectively.
    '''
   Statspec = cLoadsFacs(Statload, Statfac)
   statcoretest = []
   statshelltest = []
   statcore = 1000
   statshell = 1000
        
   for i in range(components):
       statcoretest.append(SpecErrAbs2D(Statspec.inav[i], core))
       statshelltest.append(SpecErrAbs2D(Statspec.inav[i], shell))
       if statcoretest[i] < statcore:
        statcore = statcoretest[i]
        index_c = i
        if statshelltest[i] < statshell:
            statshell = statshelltest[i]
            index_s = i
       if statcore == statshell:
           print('Warning! It guessed the same component for core and shell.')
       return index_c,index_s #index för den som bäst passar core respektive shell

class postStatprocess:
    ''' Performs all the preprocessing that only requires the original\n
    core & shell images and factors and loadings from statistical analysis.\n
    Standard amount of components is 2 and no calibration spectrum has to be choosen but it is adviced.\n
    Creates an object with the attributes:\n
    core, shell, full = reconstructed images of phases.\n
    errtot, errcore, errshell = som of errors of each phase.\n
    errmaptot, errmapcore, errmapshell, relmaptot, relmapcore, relmapshell\n
    = absolute and relative error maps of the phases.\n
    For quantification use "quantify" function.
    '''
    def __init__(self,core,shell,Statfac,Statload,components=2,calSpec=[]):
        '''  '''
        self.calSpec = calSpec
        self.factors = Statfac
        self.loadings = Statload
        particle = core + shell
        Statspec = cLoadsFacs(Statload, Statfac)
        statcoretest = []
        statshelltest = []
        statcore = 1000
        statshell = 1000
        
        for i in range(components):
            statcoretest.append(SpecErrAbs2D(Statspec.inav[i], core))
            statshelltest.append(SpecErrAbs2D(Statspec.inav[i], shell))
            if statcoretest[i] < statcore:
                statcore = statcoretest[i]
                c = i
            if statshelltest[i] < statshell:
                statshell = statshelltest[i]
                s = i
        if statcore == statshell:
            print('Warnig! It guessed the same component for core and shell.')
                
        self.core = Statspec.inav[c]
        self.corecomp = c
        self.shell = Statspec.inav[s]
        self.shellcomp = s
        self.full = self.core + self.shell
        # Sum of Errors:
        self.errtot = SpecErrAbs2D(self.full,particle)
        self.errcore = SpecErrAbs2D(self.core,core)
        self.errshell = SpecErrAbs2D(self.shell,shell)
        
        # Error maps:
        errmaptot = specMapDiff(self.full, particle)
        errmapcore = specMapDiff(self.core, core)
        errmapshell = specMapDiff(self.shell, shell)
        if len(calSpec)==0:
            print('Warning! Not specifying a calSpec can mess with the error maps.')
            errmaptot = errmaptot
            errmapcore = errmapcore
            errmapshell = errmapshell
            imparticle = particle
            
        else:
            errmaptot = setCalibration(errmaptot, calSpec)
            errmapcore = setCalibration(errmapcore, calSpec)
            errmapshell = setCalibration(errmapshell, calSpec)
            imparticle = setCalibration(particle, calSpec)
        
        self.errmaptot = errmaptot.get_lines_intensity()
        self.errmapcore = errmapcore.get_lines_intensity()
        self.errmapshell = errmapshell.get_lines_intensity()
        imComp = imparticle.get_lines_intensity()
        
        self.relmaptot = [rel(self.errmaptot[0], imComp[0]),rel(self.errmaptot[1], imComp[1])]
        self.relmapcore = [rel(self.errmapcore[0], imComp[0]),rel(self.errmapcore[1], imComp[1])]
        self.relmapshell = [rel(self.errmapshell[0], imComp[0]),rel(self.errmapshell[1], imComp[1])]
        
        # Quantification
    def quantify(self,kfac):
        ''' Quantifies the elements in the core and shell using choosen k-factors. '''
        factors = self.factors
        factors = setCalibration(factors, self.calSpec)
    
        bw = factors.inav[self.corecomp].estimate_background_windows(line_width=[5.0, 7.0])
        intensities = factors.inav[self.corecomp].get_lines_intensity(background_windows=bw)
        self.coreCu = factors.inav[self.corecomp].quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data
        self.coreAg = factors.inav[self.corecomp].quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data
        bw = factors.inav[self.shellcomp].estimate_background_windows(line_width=[5.0, 7.0])
        intensities = factors.inav[self.shellcomp].get_lines_intensity(background_windows=bw)
        self.shellCu = factors.inav[self.shellcomp].quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data
        self.shellAg = factors.inav[self.shellcomp].quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data
        
        