# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 13:47:13 2021

@author: Jonas
"""
%matplotlib qt
import hyperspy.api as hs
import numpy as np
from coreshellp import CoreShellP
from edxmat import EdxMat

#x=CoreShellP(50,20.0,15.0,100**(-3),100**(-3),1)
#core=x.core;
#shell=x.shell             #Skapar kärnan (Cu)
x1=EdxMat(50,20.0,10**(-3),1.0)  #Den första hela "partikeln" som sen kommer bli skal (Ag)
x2=EdxMat(50,15.0,10**(-3),1.0)  #Den andra mindre "partikeln" som kommer ta bort insidan så det blir ett skal (Ag)
x3=EdxMat(50,15,10**(-3),1)       #Den faktiska mindre "partikeln" som är koppar kärnan (Cu)
shell = x1.mat-x2.mat         #Skapar skalet (Ag)
core = x3.mat                 #Skapar kärnan (Cu)

'''
Det var tydligen rätt viktigt att spectrumen var gjorda på samma sätt för att man skulle kunna se skal och kärna i 
samma bild. Annars blev det att kärnan eller skalet fanns där men inte syntes varken i bilden eller spektrat då
intensiteten var för låg. 
'''
sAg = hs.load("./Spectra/PureAgFilip.msa",signal_type="EDS_TEM")
sCu = hs.load("./Spectra/PureCuFilip.msa",signal_type="EDS_TEM") 

#Tomma matriser med rätt dimensioner
arr = np.empty((50,50,2048)) 
core3d = arr.copy() # .copy() för att annars pekar de på varandra och core och shell blandas
shell3d = arr.copy()
Ag = sAg.data #Extraherar bara spectrumen från .msa-filerna
Cu = sCu.data


for i in range(0,50):
    for j in range(0,50):
        #Här fylls de tomma core och shell matriserna med spectrum*intensitet från sfärmatriserna
        core3d[i,j,0:2048]=core[i,j]*Cu
        shell3d[i,j,0:2048]=shell[i,j]*Ag
        
        
#Nu slås de core och shell ihop och blir den fullständiga particeln:
matr = shell3d + core3d
maxval = np.max(matr)
#matr = (1/maxval)*matr;
p = hs.signals.Signal1D(matr)#Läser som HyperSpy-signal
p.add_poissonian_noise(keep_dtype=True);
p.data=(1/np.max(p.data))*p.data
#Ett gäng metadata som kanske inte är nödvändigt men underlättar
p.set_signal_type("EDS_TEM") 
p.axes_manager.signal_axes[0].units = 'keV' #OBS! Enheten på energiaxeln är viktig för att kunna plotta
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.add_elements(['Ag','Cu']) #Lägger in element igen tydligen förs de inte med 
p.add_lines(['Ag_La','Cu_Ka'])

p.get_calibration_from(sAg)
p.axes_manager.indices = [25,25] #Sätter "pekaren" på mitten av partikeln, för att få spektrat från den punkten i plotten
p.plot(True)           

im = p.get_lines_intensity()
hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
kfactors = [2.32, 1.58] # [Ag_La,Cu_Ka] tagna från Boken TransmissionElectronMicroscopy (finns på LUB). Cu är för Ka, Ag La+Lb

q=p.inav[25,25].isig[0.0:10000.0]
bw = q.estimate_background_windows(line_width=[5.0, 7.0])

intensities = q.get_lines_intensity(background_windows=bw)
weight_percent = q.quantification(intensities, method='CL',
                                  factors=kfactors,composition_units='weight')
print("Vikt% Ag: "+str(weight_percent[0].data[0]))
print("Vikt% Cu: "+str(weight_percent[1].data[0]))

wtAg = shell[25,25]/(shell[25,25]+core[25,25])
wtCu = core[25,25]/(shell[25,25]+core[25,25])
print("Sann Vikt% Ag: "+str(wtAg))
print("Sann Vikt% Cu: "+str(wtCu))
q.plot(True, integration_windows='auto',background_windows=bw)

#p.decomposition(algorithm='NMF',output_dimension = 3)
#factors = p.get_decomposition_factors()
#hs.plot.plot_spectra(factors,style='cascade')