%matplotlib qt
import hyperspy.api as hs
import math
import numpy as np
import matplotlib.pyplot as plt
from edxmat import EdxMat
from coreshellp import CoreShellP, CoreShellSpec
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

#Ett gäng metadata som kanske inte är nödvändigt men underlättar
#axmanageq(matr)
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

p.decomposition(algorithm='NMF',output_dimension = 3)
factors = p.get_decomposition_factors()
hs.plot.plot_spectra(factors,style='cascade')