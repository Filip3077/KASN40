# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 13:45:25 2021

@author: Filip
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 01:11:49 2021

@author: Filip
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np
from sklearn.linear_model import LinearRegression

sCal = hs.load("./Spectra/20 nm cube Zn100Au0 SDD.msa",signal_type="EDS_TEM")

#s.inav[0] = (hs.load("../Spectra/Cu100.msa",signal_type="EDS_TEM"))
#s.inav[1] = (hs.load("../Spectra/Cu80Ag20.msa", signal_type="EDS_TEM"))
#s.inav[2] = (hs.load("../Spectra/Cu60Ag40.msa", signal_type="EDS_TEM"))
#s.inav[3] = (hs.load("../Spectra/Cu40Ag60.msa", signal_type="EDS_TEM"))
#s.inav[4] = (hs.load("../Spectra/Cu20Ag80.msa", signal_type="EDS_TEM"))
#s.inav[5] = (hs.load("../Spectra/Ag100.msa", signal_type="EDS_TEM"))

s = hs.load("./Spectra/20 nm cube Zn*.msa",stack=True,signal_type="EDS_TEM")
s.set_signal_type("EDS_TEM")
#s.get_calibration_from(sCal)
s.add_elements(['Au','Zn']) 
s.add_lines(['Zn_Ka','Au_Ma'])


I = []

for spectrum in s:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    I.append(intensities[1].data[0]/intensities[0].data[0])

I2 = np.array(I[1:-1]).reshape(-1,1)
xAu = np.linspace(0,1,11)
Ca = [x/(1-x) for x in xAu[1:-1]]

DTSAmodel = LinearRegression()
DTSAmodel.fit(I2,Ca)
DTSA_k = DTSAmodel.coef_
DTSA_factors = [1,DTSA_k[0]]
DTSAmodel.score(I2,Ca)

sLC = hs.signals.Signal1D(np.zeros((11,2048)))
sLC.set_signal_type("EDS_TEM")
sLC.get_calibration_from(s.inav[0])
sLC.add_elements(['Zn','Au']) 
sLC.add_lines(['Zn_Ka','Au_Ma'])
ILC = []

xCounter = 0
for i in range(11):
    sLC.inav[i] = s.inav[0]*(1-xCounter)+s.inav[-1]*(xCounter)
    bgLC = sLC.inav[i].estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = sLC.inav[i].get_lines_intensity(background_windows=bgLC)
    ILC.append(intensities[1].data[0]/intensities[0].data[0])
    xCounter += 0.1

sLC.set_signal_type("EDS_TEM")

sLC.add_elements(['Zn','Au']) 
sLC.add_lines(['Zn_Ka','Au_Ma'])


ILC2 = np.array(ILC[1:-1]).reshape(-1,1)
LCmodel = LinearRegression()
LCmodel.fit(ILC2,Ca)
LC_k = LCmodel.coef_
LC_factors = [1,LC_k]



DTSA_DTSA_xAg = []

for spectrum in s:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    DTSA_DTSA_xAg.append(spectrum.quantification(intensities, method='CL', factors=DTSA_factors,composition_units='weight')[1].data[0])
    
    
LC_DTSA_xAg = []

for spectrum in s:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    LC_DTSA_xAg.append(spectrum.quantification(intensities, method='CL', factors=LC_factors,composition_units='weight')[1].data[0])
    

DTSA_LC_xAg = []

for spectrum in sLC:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    DTSA_LC_xAg.append(spectrum.quantification(intensities, method='CL', factors=DTSA_factors,composition_units='weight')[1].data[0])
    

LC_LC_xAg = []

for spectrum in sLC:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    LC_LC_xAg.append(spectrum.quantification(intensities, method='CL', factors=LC_factors,composition_units='weight')[1].data[0])
    
    
DTSA_DTSA_xAg = np.array(DTSA_DTSA_xAg)/100    
LC_DTSA_xAg = np.array(LC_DTSA_xAg)/100
DTSA_LC_xAg = np.array(DTSA_LC_xAg)/100
LC_LC_xAg = np.array(LC_LC_xAg)/100

DTSA_DTSA_diff = []
LC_DTSA_diff = []
DTSA_LC_diff = []
LC_LC_diff = []

for i in range(11):
    DTSA_DTSA_diff.append(abs(DTSA_DTSA_xAg[i] - xAu[i]))
    LC_DTSA_diff.append(abs(LC_DTSA_xAg[i] - xAu[i]))
    DTSA_LC_diff.append(abs(DTSA_LC_xAg[i] - xAu[i]))
    LC_LC_diff.append(abs(LC_LC_xAg[i] - xAu[i] ))

plt.figure()
plt.subplot(211)
plt.title('DTSA spectrum')
plt.plot(xAu,DTSA_DTSA_xAg,'r',xAu,xAu,'k',xAu,LC_DTSA_xAg,'b')
plt.legend(['DTSA-k','True','Linear Combination-k'])

plt.subplot(212)
plt.plot(xAu,DTSA_DTSA_diff,'ro',xAu,LC_DTSA_diff,'bo')

plt.figure()
plt.subplot(211)
plt.title('Linear combinations')
plt.plot(xAu,DTSA_LC_xAg,'r',xAu,xAu,'k',xAu,LC_LC_xAg,'b')
plt.legend(['DTSA-k','Intended','Linear Combination-k'])

plt.subplot(212)
plt.plot(xAu,DTSA_LC_diff,'ro',xAu,LC_LC_diff,'bo')

plt.figure()
plt.plot(xAu[1:-1]/(1-xAu[1:-1]),I[1:-1],'r',xAu[1:-1]/(1-xAu[1:-1]),ILC[1:-1],'b')
