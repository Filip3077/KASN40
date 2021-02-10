# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:20:11 2021

@author: Jonas
"""

from edxmat import EdxMat

class CoreShellP:
    def __init__(self,size,r1:float,r2:float,dens1:float,dens2:float,l:float):
        '''Constructs an object modeling a core-shell spherical particle.\n
        The matrices core and shell model the core and shell respectively.\n
        All argument are defined as for EdxMat'''
        x=EdxMat(size,r1,dens1,l);
        y=EdxMat(size,r2,dens2,l);
        if r1>=r2:
            z=EdxMat(size,r2,dens1,l)
            self.shell=x.mat-z.mat;
            self.core=y.mat;
        else:
            z=EdxMat(size,r1,dens2,l)
            self.shell=y.mat-z.mat;
            self.core=x.mat;
    