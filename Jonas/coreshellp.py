# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:20:11 2021

@author: Jonas
"""
import math
import hyperspy.api as hs
import numpy as np
import scipy.misc

def cut_spectrum_bottom(hs,at):
    at=float(at);
    hs.isig[0:at].data-=hs.isig[0:at].data;
    return None

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
                    mat[i][j]=2*dens*self.__thick(i,j);
        
    def __thick(self,n,m):
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
         self.signal=signal
         L=len(spec1.data)
         arr=np.empty((a.size,a.size,L))
         core=arr.copy()
         shell = arr.copy()
         for i in range(0,a.size):
             for j in range(0,a.size):
                 core[i,j,0:L]=a.core[i,j]*spec1.data
                 shell[i,j,0:L]=a.shell[i,j]*spec2.data
         if signal:
             self.core = hs.signals.Signal1D(core)
             self.shell = hs.signals.Signal1D(shell)
         else:
             self.core=core;
             self.shell=shell;
         self.full=self.core+self.shell   
         self.base=a;
         self
     def getmatr(self):
         '''Returns the matrix representing the entire core-shell particle'''
         return self.core+self.shell
     def is_signal(self):
         return self.signal;
     def add_background(self,backspec, thickness):
         for i in range(self.base.size):
             for j in range(self.base.size):
                 if self.signal:
                     self.core.data[i][j]+=backspec.data*thickness;
                 else:
                     self.core[i][j]+=backspec.data*thickness;
                     
     def add_core_component(self, compspec,frac):
         '''Adds an element to the core. Only valid for small additions (frac<<1)'''
         for i in range(self.base.size):
            for j in range(self.base.size):
                if self.signal:
                     self.core.data[i][j]=self.core.data[i][j]*(1-frac)+compspec.data*frac*self.base.core[i,j];
            else:
                     self.core[i][j]=self.core[i][j]*(1-frac)+compspec.data*frac*self.base.core[i,j];
                     
     def add_shell_component(self,compspec,frac):
         '''Adds an element to the shell. Only valid for small additions (frac<<1)'''
         for i in range(self.base.size):
             for j in range(self.base.size):
                 if self.signal:
                     self.shell.data[i][j]=self.shell.data[i][j]*(1-frac)+compspec.data*frac*self.base.shell[i,j];
                 else:
                     self.shell[i][j]=self.shell[i][j]*(1-frac)+compspec.data*frac*self.base.shell[i,j];

    
class CoreShellBack:
    def __init__(self,a,spec,dens,signal=True):
        self.signal=signal;
        spec=spec.dens
        l=a.size()
        L=len(spec)
        mat=np.zeros((l,l))
        if signal:
            for i in range(l):
                for j in range(l):
                    mat[i][j][0:L]=hs.signals.Signal1D(spec);
        else:
            for i in range(l):
                for j in range(l):
                    mat[i][j][0:L]=spec;            
        self.size=a.size;
        self.base=a.base;
        self.back=mat
        
    def get_backg(self):
        return self.back;
    def is_signal(self):
        return self.signal