import matplotlib as mpl
import hyperspy.api as hs
import numpy as np
%matplotlib inline 
#Hade problem med att visa plots utanför själva browsern

#Ett försök att skapa en motsvarande EDS karta. I detta fall 6x6 pixlar med en energidimesion baserat på 1 genererat EDX-spektra, detta ger storleken på matrisen 6,6,2048. Har inte riktigt fått till någon skalning av spektrat i varje pixel, men det går att läsa med Hyperspys funktioner. 


k = hs.signals.Signal1D(np.random.random((6,6,2048))) #Skapar en 
s = hs.load("testspectraAu.msa",reader="msa",signal_type = "EDS_TEM")
s.set_elements(['Au'])
s.add_elements(['Cu'])
s.add_lines()
p = s+k
p.axes_manager[1].name = 'x'
p.axes_manager[0].name = 'y'

for spectra in p:
    spectra = spectra

im = p.get_lines_intensity()
hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
p.plot()
p.get_lines_intensity(['Au_Ma'],plot_result=True)
