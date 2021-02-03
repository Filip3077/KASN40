%matplotlib qt
import hyperspy.api as hs
import math
import numpy as np
import matplotlib.pyplot as plt
from edxmat import EdxMat
import scipy.misc

x1=EdxMat(50,20.0,10.49,1.0)  #Den första hela "partikeln" som sen kommer bli skal (Ag)
x2=EdxMat(50,15.0,10.49,1.0)  #Den andra mindre "partikeln" som kommer ta bort insidan så det blir ett skal (Ag)
x3=EdxMat(50,15,8.96,1)       #Den faktiska mindre "partikeln" som är koppar kärnan (Cu)
shell = x1.mat-x2.mat         #Skapar skalet (Ag)
core = x3.mat                 #Skapar kärnan (Cu)


'''
Det var tydligen rätt viktigt att spectrumen var gjorda på samma sätt för att man skulle kunna se skal och kärna i 
samma bild. Annars blev det att kärnan eller skalet fanns där men inte syntes varken i bilden eller spektrat då
intensiteten var för låg. 
'''
sAg = hs.load("PureAgFilip.msa",signal_type="EDS_TEM") #Ändrade till ett annat par av spektrum som hade samma tjocklek i DSTA-II
sCu = hs.load("PureCuFilip.msa",signal_type="EDS_TEM") 

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
#matr = (matr/maxval*255).astype(np.uint8)  #Denna behövs nog inte till hyperspy men behövs om man ska göra en egentlig bild tror jag.

p = hs.signals.Signal1D(matr) #Läser in matrisen som en hyperspysignal

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

im = p.get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
kfactors = [2.32, 1.58] # [Ag_La,Cu_Ka] tagna från Boken TransmissionElectronMicroscopy (finns på LUB). Cu är för Ka, Ag La+Lb

q=p.inav[25,25].isig[0.0:10000.0] #Tar ut ett spektra vid pixel 25,25 och skär ner på energiaxel till 0-10000eV
bw = q.estimate_background_windows(line_width=[5.0, 7.0]) #Sätter vilket fönster som bakgrunden ska tas mellan relativt till varje xray-line som är vald

intensities = q.get_lines_intensity(background_windows=bw)  #Tar fram intensiter för dessa xray lines
weight_percent = q.quantification(intensities, method='CL', #Kvantifiering med Cliff-Lorimer tar ut i viktsprocent
                                  factors=kfactors,composition_units='weight')
print("Vikt% Ag: "+str(weight_percent[0].data[0])) 
print("Vikt% Cu: "+str(weight_percent[1].data[0]))

wtAg = shell[25,25]/(shell[25,25]+core[25,25]) #Beräknar sammansättningen utfrån orginalmatris
wtCu = core[25,25]/(shell[25,25]+core[25,25])
print("Sann Vikt% Ag: "+str(wtAg)) 
print("Sann Vikt% Cu: "+str(wtCu)) 
q.plot(True, integration_windows='auto',background_windows=bw) #Plottar spektrat som används



#En snabb titt på NMF av detta, kunde snabbt se att endast två faktorer var relevanta med tidigare test. 
p.decomposition(algorithm='NMF',output_dimension = 2)
factors = p.get_decomposition_factors() #Tar ut faktorerna dvs spektrum
loadings =p.get_decomposition_loadings() #loadings är det återskapade bilderna baserat på faktorerna 

#Plottar dem båda åter igen drar ner energiaxel till max 10keV för tydlighet
hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade') 
hs.plot.plot_images(loadings, cmap='mpl_colors',
            axes_decor='off', per_row=1,
            label=['Cu Core', 'Ag Shell'],
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})

#Testar att kvantifiera core-faktorn återigen endast upp till 10keV
core = factors.inav[0].isig[0.0:10000.0]
core.add_elements(['Ag','Cu']) #Lägger in element igen tydligen förs de inte med 
core.add_lines(['Ag_La','Cu_Ka'])
corebw = core.estimate_background_windows(line_width=[5.0, 7.0]) 

intensities = core.get_lines_intensity(background_windows=corebw)
weight_percent = core.quantification(intensities, method='CL',
                                  factors=kfactors,composition_units='weight')
print("Vikt% Ag: "+str(weight_percent[0].data[0]))
print("Vikt% Cu: "+str(weight_percent[1].data[0])) 


#Väldigt likt det som fanns tidigare, vid den första kvantifieringen av grundspektrumet vid [25,25], blir nog lite att klura på om jag gjort rätt med allt
