# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 17:52:49 2021

@author: Jonas
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np
import specerr as spe
from sklearn.linear_model import LinearRegression

def silver_k_factor(hz,el,lines):
    refspec=hs.load("./Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
    s=hs.signals.Signal1D(np.zeros((11,len(hz.isig))))
    s.set_signal_type('EDS_TEM')
    s.get_calibration_from(refspec)
    
    ael=['Ag']
    for st in el:
        ael.append(st)
        
        alines=['Ag_La']
        for st in lines:
            alines.append(st)
            s.add_elements(ael)
            s.add_lines(alines)
            for i in range(len(s.isig[0])):
                s.inav[i]=hz*i*0.1+refspec*(10-i)*0.1
                #q.set_signal_type("EDS_TEM") 
                #q.add_elements(ael) 
                #q.add_lines(alines)  
    
    
    I = []
    for i in range(len(s.isig[0])):
        corebw = s.inav[i].estimate_background_windows(line_width=[2.0, 8.0]) 
        intensities = s.inav[i].get_lines_intensity(background_windows=corebw)
        I.append(intensities[1].data[0]/intensities[0].data[0])
#%%Uppskatta k-faktor
    I = I[1:len(s.isig[0])-1]
    I = np.asarray(I)
    I = I.reshape((-1,1))
    comp = np.linspace(0.1,0.9,9)
    c = comp/(1-comp)
    
    model = LinearRegression()
    model.fit(I,c)
    k = model.coef_
    return k