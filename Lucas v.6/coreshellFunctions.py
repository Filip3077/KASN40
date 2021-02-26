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
        factors = self.factors.data
        co = hs.signals.Signal1D(factors[self.corecomp])
        # co = setCalibration(co,self.calSpec)
        sh = hs.signals.Signal1D(factors[self.shellcomp])
        # sh = setCalibration(sh,self.calSpec)
        
        co.set_signal_type("EDS_TEM")
        sh.set_signal_type("EDS_TEM")
        co.axes_manager.signal_axes[0].units = 'keV'
        sh.axes_manager.signal_axes[0].units = 'keV'
        co.add_elements(['Ag','Cu'])
        sh.add_elements(['Ag','Cu'])
        co.add_lines(['Ag_La','Cu_Ka'])
        sh.add_lines(['Ag_La','Cu_Ka'])
        co.get_calibration_from(self.calSpec)
        sh.get_calibration_from(self.calSpec)
        bw = co.estimate_background_windows(line_width=[5.0, 2.0])
        intensities = co.get_lines_intensity(background_windows=bw)
        self.coreCu = co.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data
        self.coreAg = co.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data
        intensities = sh.get_lines_intensity(background_windows=bw)
        self.shellCu = sh.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data
        self.shellAg = sh.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data
        
        
    
    
    
    
    
        
    


