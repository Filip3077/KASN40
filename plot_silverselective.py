import matplotlib as mpl

import hyperspy.api as hs
import numpy as np
%matplotlib qt 
#Hade problem med att visa plots utanför själva browsern

#Ett försök att skapa en motsvarande EDS karta. I detta fall 6x6 pixlar med en energidimesion baserat på 1 genererat EDX-spektra, detta ger storleken på matrisen 6,6,2048. Har inte riktigt fått till någon skalning av spektrat i varje pixel, men det går att läsa med Hyperspys funktioner. 

k = np.ones((6,6,2048)) #Skapar grundmatrisen 
s = hs.load("PureAg.msa",reader="msa",signal_type = "EDS_TEM") #Läser in spektra med rent silver (skapat i DSTA-II)


s.add_elements(['Ag']) #Lägger till element till spektrat (kanske kan tas bort)
s.add_lines() #Motsvarande x-ray lines

#Loopar igenom matrisen k och för varje pixel i 6x6 planet skriver jag in spektrumet med en skalfaktor baserat på x
cols = 6
rows = 6
for x in range(rows):
    for y in range(cols):
        k[x,y,:] = s.data*x
        

p = hs.signals.Signal1D(k) #Läser in matrisen som en hyperspysignal

#Ett gäng metadata som kanske inte är nödvändigt men underlättar
p.set_signal_type("EDS_TEM") 
p.axes_manager.signal_axes[0].units = 'keV' #OBS! Enheten på energiaxeln är viktig för att kunna plotta
p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'
p.add_elements(['Ag']) #Lägger in element igen tydligen förs de inte med 
p.add_lines(['Ag_La'])

im = p.get_lines_intensity() #Här kommer problemen går att plotta vanligt och detta ger ett resultat som man kan förvänta sig. Men här kommer probelm med DataAxis något med indexet på energiaxeln tror jag. 
#get_line_intesity ska ta de inlagda xray-lines och endast plotta dessa i varje pixel, dessa verkar hittas utan det är något i datastrukturen hos p som ej funkar 


#Kod för att plotta en färgkodad bild direkt från dokumentationen
#hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
  #  colorbar='single', vmin='1th', vmax='99th', scalebar='all',
 #   scalebar_color='black', suptitle_fontsize=16,
 #   padding={'top':0.8, 'bottom':0.10, 'left':0.05,
 #            'right':0.85, 'wspace':0.20, 'hspace':0.10})
