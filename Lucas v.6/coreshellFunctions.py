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

class postStatprocess:
    '''  '''
    def __init__(self,core,shell,Statfac=[],Statload=[],calSpec=[]):
        '''  '''
        particle = core + shell
        Statspec = cLoadsFacs(Statload, Statfac)
        Statcoretest1 = SpecErrAbs2D(Statspec.inav[0],core)
        Statcoretest2 = SpecErrAbs2D(Statspec.inav[1],core)
        Statshelltest1 = SpecErrAbs2D(Statspec.inav[1],shell)
        Statshelltest2 = SpecErrAbs2D(Statspec.inav[0],shell)
        print('Coretest1: '+str(Statcoretest1)+'\nCoretest2: '+str(Statcoretest2)+'\nShelltest1: '+str(Statshelltest1)+'\nShelltest2: '+str(Statshelltest2))
        if Statcoretest1 < Statcoretest2 and Statshelltest1 < Statshelltest2:
            self.core = Statspec.inav[0]
            self.shell = Statspec.inav[1]
            self.full = self.core + self.shell
        elif Statcoretest1 > Statcoretest2 and Statshelltest1 > Statshelltest2:
            self.core = Statspec.inav[1]
            self.shell = Statspec.inav[0]
            self.full = self.core + self.shell
        else:
            print('One or more of the reconstructed images match neither the core nor shell')
        
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
        
        
        
    
    
    
    
    
        
    


