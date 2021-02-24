# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:07:34 2021

@author: Lucas
"""
import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec
from specMapDiff import *
from specerr import *

def genfullparticle(size,r1,r2,specCore,specShell,signal=True,dens1=1,dens2=1,l=1):
    ''' Generates a full core@shell particle using coreshellp.py \n
    Only size,r1,r2 and the spectrums is necessary, densities and pixel size \n
    is set to 1 by defult but can be changed.\n
    Whichever of r1 or r2 can be the large or small radius.'''
    a = CoreShellP(size,r1,r2,dens1,dens2,l) 
    return CoreShellSpec(a,specCore,specShell,signal)

class
def postNMFBSSprocess(core,shell,NMFfac=[],NMFload=[],BSSfac=[],BSSload=[]):
    '''  '''
    if len(NMFfac) == 0:
        
    else:
        NMFspec = cLoadsFacs(loadings, factors)
        NMFcoretest1 = SpecErrAbs2D(NMFspec.inav[0],core)
        NMFcoretest2 = SpecErrAbs2D(NMFspec.inav[1],core)
        if NMFcoretest1 > NMFcoretest2:
            NMFcore = NMFspec.inav[0]
            NMFshell = NMFspec.inav[1]
        else:
            NMFcore = NMFspec.inav[1]
            NMFshell = NMFspec.inav[0]
    if len(BSSfac) == 0:
        
    else:
        BSSspec = cLoadsFacs(bssload,bssfac)
        BSScoretest1 = SpecErrAbs2D(BSS.inav[0],core)
        BSScoretest2 = SpecErrAbs2D(BSS.inav[1],core)
        if BSScoretest1 > BSScoretest2:
            BSScore = BSSspec.inav[0]
            BSSshell = BSSspec.inav[1]
        else:
            BSScore = BSSspec.ina[1]
            BSSshell = BSSspec.inav[0]
        
    
    
    
    
    
        
    


