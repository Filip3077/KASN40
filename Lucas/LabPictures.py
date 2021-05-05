# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 09:03:45 2021

@author: Lucas
"""

import hyperspy.api as hs
from coreshellfunctions import*

#%%

s1 = hs.load('./images/EDS Data_particle 1.rpl').as_signal1D(0)
s2 = hs.load('./images/EDS Data_particle 2.rpl').T
s3 = hs.load('./images/EDS Data_particle 3.rpl').T
cal = hs.load('../Spectra/20nm cube Cu100Ag0.msa')

s1 = AuZnCalibration(s1,1/6,['Au','Zn','Fe'])
s2 = AuZnCalibration(s2)
s3 = AuZnCalibration(s3)

s1 = s1.rebin(scale=[4,4,1])
s2 = s2.rebin(scale=[4,4,1])
s3 = s3.rebin(scale=[4,4,1])

#%%

redbluePlot(s1,'AuZn_22nm 1'); redbluePlot(s2,'AuZn_22nm 2'); redbluePlot(s3,'AuZn_22nm 3')

#%%

s1.decomposition(True); s2.decomposition(True); s3.decomposition(True)

s1.plot_explained_variance_ratio(); s2.plot_explained_variance_ratio(); s3.plot_explained_variance_ratio()
nfac = 3

s1.decomposition(True, algorithm='NMF', output_dimension=nfac)
s2.decomposition(True, algorithm='NMF', output_dimension=nfac)
s3.decomposition(True, algorithm='NMF', output_dimension=nfac)

factors = []; loadings = []
factors.append(s1.get_decomposition_factors().inav[0:nfac])
factors.append(s2.get_decomposition_factors().inav[0:nfac])
factors.append(s3.get_decomposition_factors().inav[0:nfac])
loadings.append(s1.get_decomposition_loadings().inav[0:nfac])
loadings.append(s2.get_decomposition_loadings().inav[0:nfac])
loadings.append(s3.get_decomposition_loadings().inav[0:nfac])

#%%
facs = hs.stack(factors[0:nfac+1])
loads = [loadings[0],loadings[1],loadings[2]]
hs.plot.plot_spectra(facs,style='cascade')
plt.title('PCA')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)    
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    
hs.plot.plot_images(loads, cmap='mpl_colors',
                    axes_decor='off', per_row=nfac,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%

def quantify(factors, kfac):
    ''' Quantifies the elements in the core or shell using choosen k-factors. '''
    # factors = setCalibration(factors, calSpec)

    bw = factors.estimate_background_windows(line_width=[5.0, 7.0])
    intensities = factors.get_lines_intensity(background_windows=bw)
    wtA = factors.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data
    wtB = factors.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data
    return wtA, wtB

kfac = [2,585,1,401]
# wtAu1, wtZn1 = quantify(s1.inav[:,:].sum(),kfac)

bw = s1.inav[:,:].sum().estimate_background_windows(line_width=[3.0, 2.0])

s1.inav[:,:].sum().plot(background_windows=bw)

# s.get_lines_intensity(background_windows=bw, plot_result=True)


