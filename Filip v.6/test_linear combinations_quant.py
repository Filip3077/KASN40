# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:53:32 2021

@author: Jonas Elmroth Nordlander
KASN20 Projekt VT 2021
Skriptet testar hur väl det funkar att linjärkombinera spektrum för 100 %\n
rent Cu och Ag (både i molprocent och mass procent) jämtemot MC-simulerade \n
spektrum med motsvarande sammansättning.

"""
import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np
import specerr as spe
from sklearn.linear_model import LinearRegression



#Importera spektra av en kub med angiven sammansättning.
MCu=63.55;#Mw Cu g/mol
MAg=107.9;#Mw Ag g/mol
s=np.zeros((1,11,2048));
sAg = hs.load("../Spectra/MCSim 100% Cu.msa",signal_type="EDS_TEM")

s[0][0]=hs.load("../Spectra/MCSim 100% Cu.msa",signal_type="EDS_TEM")
s[0][1]=hs.load("../Spectra/MCSim 90Cu10Ag.msa", signal_type="EDS_TEM")
s[0][2]=hs.load("../Spectra/MCSim 80Cu20Ag.msa", signal_type="EDS_TEM");
s[0][3]=hs.load("../Spectra/MCSim 70Cu30Ag.msa", signal_type="EDS_TEM");
s[0][4]=hs.load("../Spectra/MCSim 60Cu40Ag.msa", signal_type="EDS_TEM");
s[0][5]=hs.load("../Spectra/MCSim 50Cu50Ag.msa", signal_type="EDS_TEM");
s[0][6]=hs.load("../Spectra/MCSim 40Cu60Ag.msa", signal_type="EDS_TEM");
s[0][7]=hs.load("../Spectra/MCSim 30Cu70Ag.msa", signal_type="EDS_TEM");
s[0][8]=hs.load("../Spectra/MCSim 20Cu80Ag.msa", signal_type="EDS_TEM");
s[0][9]=hs.load("../Spectra/MCSim 10Cu90Ag.msa", signal_type="EDS_TEM");
s[0][10]=hs.load("../Spectra/MCSim 100% Ag.msa", signal_type="EDS_TEM");


#Konverterar till ett Hyperspy objekt för att kunna lägga på metadata
p = hs.signals.Signal1D(s)
p.set_signal_type("EDS_TEM") 
p.axes_manager.signal_axes[0].units = 'keV' 
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.add_elements(['Ag','Cu']) 
p.add_lines(['Ag_La','Cu_Ka'])  
p.get_calibration_from(sAg)


#Baserat på de spektrumen genererade av DSTA-II beräknas IAg/ICu för alla spektrum 
I = []
for i in p:
    corebw = i.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = i.get_lines_intensity(background_windows=corebw)
    I.append(intensities[0].data[0]/intensities[1].data[0])
    
    
#Gör en passning för k-faktor till dessa intensiteter
I = I[1:10]
I = np.asarray(I)
I = I.reshape((-1,1))
comp = np.linspace(0.1,0.9,9)
c = comp/(1-comp)

model = LinearRegression()
model.fit(I,c)
k = model.coef_
kfactors = [k, 1]

#Med denna faktor beräknas sedan At% via CL
factors1 = []
for i in p:
    corebw = i.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = i.get_lines_intensity(background_windows=corebw)
    weight_percent = p.quantification(intensities, method='CL', factors=kfactors)
    factors1.append(weight_percent[0].data[0])

abserr=np.zeros(9);
sqerr=np.zeros(9);
q = []

#En linjärkombination motsvarande de simulerade spektrumen från DSTA
for i in range(1,10):
    spec2=(1-.1*i)*s[0][0]+.1*i*s[0][-1];
    abserr[i-1]=spe.SpecErrAbs(s[0][i],spec2);
    sqerr[i-1]=spe.SpecErrSq(s[0][i],spec2);
    tempspec = hs.signals.Signal1D(spec2)
    q.append(tempspec)


AgKvant = []

I2 = []
#Tar fram intensiterna för Ag och Cu för dessa spektrum Alltså testar om det blir någon skillnad i k-faktorn om den är passad till 
#linjärkombinationen eller genererade spektra från DSTA
for i in q:
    i.set_signal_type("EDS_TEM") 
    i.axes_manager.signal_axes[0].units = 'keV'
    i.add_elements(['Ag','Cu']) 
    i.add_lines(['Ag_La','Cu_Ka'])  
    i.get_calibration_from(sAg)
    corebw = i.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = i.get_lines_intensity(background_windows=corebw)
    weight_percent = i.quantification(intensities, method='CL', factors=kfactors)
    AgKvant.append(weight_percent[0].data[0])
    I2.append(intensities[0].data[0]/intensities[1].data[0])

#Likt tidigare passas en k-faktor till denna data. 
I2 = np.asarray(I2)

I2 = I2.reshape((-1,1))

model2 = LinearRegression()
model2.fit(I2,c)
k2 = model2.coef_
kfactors2 = [k2, 1]


#Beräknar tillbaka At% silver igen 
AgKvantV2 = []
for i in q:
    i.set_signal_type("EDS_TEM") 
    i.axes_manager.signal_axes[0].units = 'keV' 
    i.add_elements(['Ag','Cu']) 
    i.add_lines(['Ag_La','Cu_Ka'])  
    i.get_calibration_from(sAg)
    corebw = i.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = i.get_lines_intensity(background_windows=corebw)
    weight_percent = i.quantification(intensities, method='CL', factors=kfactors2)
    AgKvantV2.append(weight_percent[0].data[0])

#Beräknar At% silver för DTSA spektrumen med faktorerna från linjärkombinationbaserade faktorerna 
factors2 = []
for i in p:
    corebw = i.estimate_background_windows(line_width=[5.0, 7.0]) 
    intensities = i.get_lines_intensity(background_windows=corebw)
    weight_percent = p.quantification(intensities, method='CL', factors=kfactors2)
    factors2.append(weight_percent[0].data[0])


#Plottar alltsammans 
comp = comp*100
diff = abs(AgKvant-comp)
diff2 = abs(AgKvantV2-comp)
plt.figure()

plt.subplot(211)
plt.title('Measured At% vs True, on linear combination')
plt.plot(comp,AgKvant,'b', comp,comp,'k',comp,AgKvantV2,'r')
plt.legend(["Facs based on DSTA","True values","Facs based on linear comb."])
plt.xlabel('At% Ag true')
plt.ylabel('At% Ag calculated')

plt.subplot(212)
plt.title('Diff between true and calculated')
plt.plot(comp,diff,'bo',comp,diff2,'ro')
plt.xlabel('At% Ag')

comp = np.linspace(0,1,11)
factors1 = [x/100 for x in factors1]
factors2 = [x/100 for x in factors2]

plt.figure()
plt.title('Measured At% vs True, on DSTA spectra')
plt.legend(["Facs based on DSTA","Facs based on linear comb.","True values"])
plt.plot(comp,factors1,'b',comp,factors2,'r',comp,comp,'k')
plt.xlabel('At% Ag')
