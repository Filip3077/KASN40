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

sCal = hs.load("../Spectra/Cu100.msa",signal_type="EDS_TEM")

#s.inav[0] = (hs.load("../Spectra/Cu100.msa",signal_type="EDS_TEM"))
#s.inav[1] = (hs.load("../Spectra/Cu80Ag20.msa", signal_type="EDS_TEM"))
#s.inav[2] = (hs.load("../Spectra/Cu60Ag40.msa", signal_type="EDS_TEM"))
#s.inav[3] = (hs.load("../Spectra/Cu40Ag60.msa", signal_type="EDS_TEM"))
#s.inav[4] = (hs.load("../Spectra/Cu20Ag80.msa", signal_type="EDS_TEM"))
#s.inav[5] = (hs.load("../Spectra/Ag100.msa", signal_type="EDS_TEM"))

s = hs.load("../Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")
s.set_signal_type("EDS_TEM")
#s.get_calibration_from(sCal)
s.add_elements(['Ag','Cu']) 
s.add_lines(['Ag_La','Cu_Ka'])


I = []

for spectrum in s:
    bg = spectrum.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = spectrum.get_lines_intensity(background_windows=bg)
    I.append(intensities[1].data[0]/intensities[0].data[0])

I2 = np.array(I[1:-1]).reshape(-1,1)
xAg = np.linspace(0,1,6)
Ca = [x/(1-x) for x in xAg[1:-1]]

DTSAmodel = LinearRegression()
DTSAmodel.fit(I2,Ca)
DTSA_k = DTSAmodel.coef_
DTSA_factors = [1,DTSA_k]
DTSAmodel.score(I2,Ca)

sLC = hs.signals.Signal1D(np.zeros((6,2048)))
sLC.set_signal_type("EDS_TEM")
sLC.get_calibration_from(s.inav[0])
sLC.add_elements(['Ag','Cu']) 
sLC.add_lines(['Ag_La','Cu_Ka'])
ILC = []

xCounter = 0
for i in range(6):
    sLC.inav[i] = s.inav[0]*(1-xCounter)+s.inav[5]*(xCounter)
    bgLC = sLC.inav[i].estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = sLC.inav[i].get_lines_intensity(background_windows=bgLC)
    ILC.append(intensities[1].data[0]/intensities[0].data[0])
    xCounter += 0.2

sLC.set_signal_type("EDS_TEM")

sLC.add_elements(['Ag','Cu']) 
sLC.add_lines(['Ag_La','Cu_Ka'])


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

for i in range(6):
    DTSA_DTSA_diff.append(abs(DTSA_DTSA_xAg[i] - xAg[i]))
    LC_DTSA_diff.append(abs(LC_DTSA_xAg[i] - xAg[i]))
    DTSA_LC_diff.append(abs(DTSA_LC_xAg[i] - xAg[i]))
    LC_LC_diff.append(abs(LC_LC_xAg[i] - xAg[i] ))

plt.figure()
plt.subplot(211)
plt.title('DTSA spectrum')
plt.plot(xAg,DTSA_DTSA_xAg,'r',xAg,xAg,'k',xAg,LC_DTSA_xAg,'b')
plt.legend(['DTSA-k','True','Linear Combination-k'])

plt.subplot(212)
plt.plot(xAg,DTSA_DTSA_diff,'ro',xAg,LC_DTSA_diff,'bo')

plt.figure()
plt.subplot(211)
plt.title('Linear combinations')
plt.plot(xAg,DTSA_LC_xAg,'r',xAg,xAg,'k',xAg,LC_LC_xAg,'b')
plt.legend(['DTSA-k','True','Linear Combination-k'])

plt.subplot(212)
plt.plot(xAg,DTSA_LC_diff,'ro',xAg,LC_LC_diff,'bo')

plt.figure()
plt.plot(xAg[1:-1]/(1-xAg[1:-1]),I[1:-1],'r',xAg[1:-1]/(1-xAg[1:-1]),ILC[1:-1],'b')
