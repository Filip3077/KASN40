# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 19:39:12 2021

@author: Jonas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt

def redBlueMap(im,supTitle=None, label=None):
    if (supTitle==None and label==None):
        hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
                            colorbar='single', vmin='1th', vmax='99th', scalebar='all',
                            scalebar_color='black', suptitle_fontsize=16,
                            padding={'top':0.8, 'bottom':0.10, 'left':0.05,
                                     'right':0.85, 'wspace':0.20, 'hspace':0.10})
    else:
      hs.plot.plot_images(im, suptitle=supTitle, label=label,tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
                            colorbar='single', vmin='1th', vmax='99th', scalebar='all',
                            scalebar_color='black', suptitle_fontsize=16,
                            padding={'top':0.8, 'bottom':0.10, 'left':0.05,
                                     'right':0.85, 'wspace':0.20, 'hspace':0.10})  
    return None

def orBlueMapCuAg(factors,loadings, title):
    hs.plot.plot_spectra(factors.isig[0.0:10000.0],style='cascade')
    plt.title(title)
    plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
    plt.axvline(8040, c='k', ls=':', lw=0.5)
    plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
    plt.axvline(930, c='k', ls=':', lw=0.5)
    plt.axvline(2984, c='k', ls=':', lw=0.5)
    plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    hs.plot.plot_images(loadings, cmap='mpl_colors',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})
    return None

def silver_k_factor(hz,el,lines):
    refspec=hs.load("20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
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