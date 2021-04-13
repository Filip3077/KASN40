# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:56:44 2021

@author: Jonas
"""

%matplotlib qt
import hyperspy.api as hs
import math
import numpy as np
import matplotlib.pyplot as plt
from edxmat import EdxMat
from coreshellp import CoreShellP, CoreShellSpec
import scipy.misc


x=CoreShellP(50,20.0,15.0,10.49,8.96,1)#Skapar ett objekt med core- och shellmatriser

sAg = hs.load("./Spectra/PureAgFilip.msa",signal_type="EDS_TEM")
sCu = hs.load("./Spectra/PureCuFilip.msa",signal_type="EDS_TEM") 
Ag = sAg.data #Extraherar bara spectrumen från .msa-filerna
Cu = sCu.data
matr=CoreShellSpec(x,Cu,Ag).getmatr()