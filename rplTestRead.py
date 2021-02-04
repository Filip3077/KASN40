%matplotlib qt
import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np

s = hs.load("CuAg850_SOI4.rpl",reader="rpl",signal_type = "EDS_TEM").T

#Tar ut den viktiga delen av bilden
p = s.inav[107:185,21:161]

cal = hs.load("AgCu.msa",signal_type="EDS_TEM")

#Även Martins verkade inte ha någon kalibrering så importerade en från en DSTA-II spektrum
p.get_calibration_from(cal)

p.axes_manager[0].name = 'y'
p.axes_manager[1].name = 'x'
p.axes_manager['x'].units = 'nm'
p.axes_manager['y'].units = 'nm'
p.axes_manager[-1].name = 'E'

#Alla spektrum var förskjutna med 200eV så lade till en offset
p.axes_manager['E'].offset = -200

 #Lägger in alla element jag hittat i bilden
p.add_elements(['Ag','Cu','C','Si'])
p.add_lines(['Ag_La','Cu_Ka','C_Ka','Si_Ka'])

#Plottar dessa 
im = p.get_lines_intensity()
hs.plot.plot_images(im,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})

#Sätter antalet komponenter jag vill ha i decompositionen, tror egentligen bara 2 är relevanta men tog med en extra
decomp_dim = 3

#Var också tvungen att byta datatyp, kör sedan både först NMF och sedan BSS med fastICA
p.change_dtype('float64')
p.decomposition(algorithm='NMF',output_dimension = decomp_dim)
factors = p.get_decomposition_factors() #Tar ut faktorerna dvs spektrum
loadings =p.get_decomposition_loadings() #loadings är det återskapade bilderna baserat på faktorerna 
p.blind_source_separation(number_of_components=decomp_dim)#,algorithm="orthomax"
bssfac = p.get_bss_factors()
bssload = p.get_bss_loadings()

 #Normerar datan så det är lättare att se vad som händer i spektrumen rent kvalitativt
for f in factors:
    f.data /= f.data.max()
    
#En massa plottande, skär av energiaxeln vid 10keV verkar inte hända så mycket där ovanför
hs.plot.plot_spectra(factors.isig[:10000.0],style='cascade',padding=-1) 

plt.title('NMF')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')
plt.axvline(277, c='k', ls=':', lw=0.5)
plt.text(x=277, y=0.8, s='C-K$_\\alpha$', color='k')
plt.axvline(1740, c='k', ls=':', lw=0.5)
plt.text(x=1740, y=0.8, s='Si-K$_\\alpha$', color='k')

hs.plot.plot_images(loadings, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})


#Plottar motvarande fast nu för BSS/fastICA
hs.plot.plot_spectra(bssfac.isig[0.0:10000.0],style='cascade') 
plt.title('NMF+BSS with fastICA')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=8040, y=1.6, s='Cu-K$_\\alpha$', color='k')
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=1.6, s='Ag-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.text(x=930, y=1.6, s='Cu-L$_\\alpha$', color='k')
plt.axvline(277, c='k', ls=':', lw=0.5)
plt.text(x=277, y=1.6, s='C-K$_\\alpha$', color='k')
plt.axvline(1740, c='k', ls=':', lw=0.5)
plt.text(x=1740, y=1.6, s='Si-K$_\\alpha$', color='k')

hs.plot.plot_images(bssload, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})



avgcounts = p.inav[:,:].data.sum()/(p.data.shape[0]*p.data.shape[1])
print("Medelantal counts: "+str(avgcounts))
#################################################
#Har inte lyckats få någon separation mellan kärna och skal endast mellan bakgrund och partikel. Kommer nog krävas en del jobb och vi ska ha så här lågt antal counts 
#Det som BSS hjälpte med verkade främst vara att fullständigt separera kol och kiselbakgrunden men än så länge inget mellan koppar och silver 
