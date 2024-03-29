# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:18:03 2021

@author: Jonas
"""

import hyperspy as hp
import numpy as np

def SpecErrAbs (spec1, spec2):
    '''Compares two spectra as vectors and generates a fractional difference between them as \n
    a percentage error of spec1'''#Input är alltså tänkt att vara på formen hp.load(fil).data
    l1=len(spec1);
    l2=len(spec2);
    s=0;
    if (l1>=l2):#Om spektrumen har olika längd så klipper funktionen av det längre spektrumet
        l=l2;#så att det matchar det kortare spektrumet.
    else:
        l=l1;
    for i in range(l):
        s=s+abs(spec1[i]-spec2[i]);#Absolutbelopp så inte negativa och positiva fel tar ut varandra.
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
    return 1-s/z;#Returnerar alltså R^2 där spec2 tolkas som en modell av spec1.