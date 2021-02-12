# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:20:11 2021

@author: Jonas
"""
import math
import hyperspy.api as hs
import numpy as np
from edxmat import EdxMat
import scipy.misc

class CoreShellP:
    def __init__(self,size,r1:float,r2:float,dens1:float,dens2:float,l:float):
        '''Constructs an object modeling a core-shell spherical particle.\n
        The object properties core and shell  are matrices \n
        modeling the core and shell respectively.\n
        All arguments are defined as for EdxMat.'''
        #Klassen reproducerar arbetet med att skapa skal- och kärnmatrisen och sparar
        #båda i ett objekt där de kan anropas med self.shell och self.core.
        x=EdxMat(size,r1,dens1,l);
        y=EdxMat(size,r2,dens2,l);
        self.size=size;
        self.l=l;
        if r1>=r2:
            z=EdxMat(size,r2,dens1,l)
            self.shell=x.mat-z.mat;
            self.core=y.mat;
            # Normaliserar:
            # maxval=np.max([np.max(self.core),np.max(self.shell)])
            # self.shell=self.shell/maxval
            # self.core=self.core/maxval
        else:
            z=EdxMat(size,r1,dens2,l)
            self.shell=y.mat-z.mat;
            self.core=x.mat;
            # Normaliserar:
            # maxval=np.max([np.max(self.core),np.max(self.shell)])
            # self.shell=self.shell/maxval
            # self.core=self.core/maxval

class CoreShellSpec:
     def __init__(self,a,specCore,specShell):
         '''Takes a CoreShellP object a and inserts the core spectrum \n
         spec1 and the shell spectrum spec2 at the appropriate places.\n
         Both spectra must be HyperSpy 1D signals.'''
         #Klassen reproducerar i princip delar av det manuella arbetet i core@shell.py för
         #att skapa själva matrisen som skall bli HyperSpy-signalen.
         #För att skapa den slutgiltiga matrisen som skall in i hs.signals.Signal1D()
         #skapar man en ny matris x=self.core+self.shell.
         L=len(specCore.data)
         arr=np.empty((a.size,a.size,L))
         core=arr.copy()
         shell = arr.copy()
         for i in range(0,a.size):
             for j in range(0,a.size):
                 core[i,j,0:L]=a.core[i,j]*specCore.data
                 shell[i,j,0:L]=a.shell[i,j]*specShell.data
         self.core=core;
         self.shell=shell;