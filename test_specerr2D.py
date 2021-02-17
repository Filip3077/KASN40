# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:00:50 2021

@author: Jonas
"""
import numpy as np
from specerr import *
import hyperspy.api as hs
import math
import matplotlib.pyplot as plt
from edxmat import EdxMat
from coreshellp import *
import scipy.misc


x=CoreShellP(50,20.0,15.0,10.49,8.96,1)#Skapar ett objekt med core- och shellmatriser
core=x.core;
shell=x.shell

'''
Det var tydligen rätt viktigt att spectrumen var gjorda på samma sätt för att man skulle kunna se skal och kärna i 
samma bild. Annars blev det att kärnan eller skalet fanns där men inte syntes varken i bilden eller spektrat då
intensiteten var för låg. 
'''
sAg = hs.load("./Spectra/PureAgFilip.msa",signal_type="EDS_TEM")
sCu = hs.load("./Spectra/PureCuFilip.msa",signal_type="EDS_TEM") 

Ag = sAg.data #Extraherar bara spectrumen från .msa-filerna
Cu = sCu.data

a=CoreShellSpec(x,Cu,Ag);#Skapar 3D-matriserna a.core och a.shell (ekvivalenta med tidigare core3d och shell3d)

#Nu slås de core och shell ihop och blir den fullständiga particeln:
matr=a.core+a.shell;
maxval = np.max(matr)
#matr = (matr/maxval*255).astype(np.uint8)  #Denna behövs nog inte till hyperspy men behövs om man ska göra en egentlig bild tror jag.

p = hs.signals.Signal1D(matr) #Läser in matrisen som en hyperspysignal
x=np.ones((2,2,2));
y=0.5*x;
q=SpecErrAbs2D(x,y)
print(q)
a=p.data
q2=SpecErrAbs2D(a,a*0.5)