# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:18:03 2021

@author: Jonas
"""

import hyperspy as hp
import numpy as np

def SpecErrAbs (spec1, spec2):
    '''Compares two spectra as vectors and generates a fractional difference between them as \n
    a percentage error of spec1'''
    l1=len(spec1);
    l2=len(spec2);
    s=0;
    if (l1>=l2):# If the spectra have different lengths the longest is cut to match the shorter one.  
        l=l2;
    else:
        l=l1;
    for i in range(l):
        s=s+abs(spec1[i]-spec2[i]);
    if (sum(spec1)==0):
        return 0
    else:
        return s/sum(spec1);

def SpecErrSq (spec1, spec2):
    ''' Takes two spectra as vectors and returns the coefficient of determination using spec2 as a model of spec1'''
    l1=len(spec1);
    l2=len(spec2);
    s=0;
    z=0;
    if (l1>=l2):
           l=l2;
    else:
           l=l1;
    for i in range(l):
            s=s+((spec1[i]-spec2[i]))**2;
            z=z+((spec1[i]-np.mean(spec1)))**2            
    return 1-s/z;# Returns R^2 where spec2 is seen as a model of spec1.

def SpecErrAbs2D(map1,relMap):
    '''Takes a numpy array with spectra stored in each element and generates\n
    a fractional difference between the two arrays based in the spectra. The norm \n
    is spec1.'''
    reldif = abs(map1-relMap).data.sum()/relMap.data.sum()
    return reldif
    
    # l1=spec1.shape;
    # l2=spec2.shape;
    # k0=min(l1[0],l2[0]);
    # k1=min(l1[1],l2[1]);
    # store=np.empty((k0,k1))
    # for i in range(k0):
    #     for j in range(k1):
    #         store[i][j]=SpecErrAbs(spec1[i][j],spec2[i][j])*np.sum(spec1[i][j])
    # return np.sum(store)/np.sum(spec1)
