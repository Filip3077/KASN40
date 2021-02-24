# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:20:11 2021

@author: Jonas
"""
import math
import hyperspy.api as hs
import numpy as np
import scipy.misc

class EdxMat:
    def __init__(self,size,r: float,dens:float,l: float):
        '''size: defines the 'picture' size as size x size pixels\n
        r: defines the radius of the spherical particle\n
        dens: density of the material the particle is made of\n
        l: pixel size in the same unit as r (i.e. the area is l^2)'''
        self.size=size;
        self.r=r;
        self.dens=dens;
        self.l=l;
        rred=r/l;
        self.mid= float((size-1)/2);#-1 because indices start from zero
        #The constructor now constructs a matrix using the input arguments
        mat=np.zeros((size,size));
        self.mat=mat;
        for i in range(size):
            for j in range(size):
                if rred**2>=((i-self.mid)**2+(j-self.mid)**2):
                    mat[i][j]=2*dens*self.thick(i,j);
        
    def thick(self,n,m):
        nr=n-self.mid;
        mr=m-self.mid;
        l=self.l;
        r=self.r;
        return math.sqrt(r**2-(nr)**2-(mr)**2)*l**2; 
    

class CoreShellP:
    def __init__(self,size,r1:float,r2:float,dens1:float,dens2:float,l:float):
        '''Constructs an object modeling a core-shell spherical particle.\n
        The object properties core and shell  are matrices \n
        modeling the core and shell respectively.\n
        All arguments are defined as for EdxMat.'''
        
        x=EdxMat(size,r1,dens1,l);
        y=EdxMat(size,r2,dens2,l);
        self.size=size;
        self.l=l;
        if r1>=r2:
            z=EdxMat(size,r2,dens1,l)
            self.shell=x.mat-z.mat;
            self.core=y.mat;

        else:
            z=EdxMat(size,r1,dens2,l)
            self.shell=y.mat-z.mat;
            self.core=x.mat;


class CoreShellSpec:
     def __init__(self,a,spec1,spec2,signal=True):
         '''Takes a CoreShellP object a and inserts the core spectrum \n
         spec1 and the shell spectrum spec2 at the appropriate places.\n
         Both spectra must be HyperSpy 1D signals.'''
         
         L=len(spec1.data)
         arr=np.empty((a.size,a.size,L))
         core=arr.copy()
         shell = arr.copy()
         for i in range(0,a.size):
             for j in range(0,a.size):
                 core[i,j,0:L]=a.core[i,j]*spec1.data
                 shell[i,j,0:L]=a.shell[i,j]*spec2.data
         full = core + shell
         if signal:
             self.core = hs.signals.Signal1D(core)
             self.shell = hs.signals.Signal1D(shell)
             self.full = hs.signals.Signal1D(full)
         else:
             self.core=core;
             self.shell=shell;
             self.full=full
         self.base=a;
     def getmatr(self):
         '''Returns the matrix representing the entire core-shell particle'''
         return self.core+self.shell 