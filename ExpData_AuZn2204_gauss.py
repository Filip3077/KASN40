# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 08:17:19 2021

@author: Filip
"""

import hyperspy.api as hs
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter;

#s20nm =  hs.load("EDS Data_particle*.rlp",stack=True,signal_type="EDS_TEM")

s1 = hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 4.rpl',signal_type="EDS_TEM").T
s2 = hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 2.rpl',signal_type="EDS_TEM").T
s3 = hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 3.rpl',signal_type="EDS_TEM").T
sRef =  hs.load("./Spectra/PureCuFilip.msa",signal_type="EDS_TEM")
sRef = sRef.isig[0:1500]
s=[s1,s2,s3];
#%%

for z in s:
    z.axes_manager[0].name = 'y'
    z.axes_manager[1].name = 'x'
    z.axes_manager[-1].name = 'E'
    z.axes_manager.signal_axes[0].units ='eV'
    z.axes_manager['E'].scale = 10
    z.axes_manager['x'].units = 'nm'
    z.axes_manager['y'].units = 'nm'
    z.axes_manager['E'].offset = -200
    z.axes_manager['y'].scale = 0.22
    z.axes_manager['x'].scale = 0.22
    z.add_elements(['Au','Zn','C','C']) #Lägger in element igen tydligen förs de inte med 
    z.add_lines(['Au_La','Zn_La','C_Ka','C_Ka'])
    z.change_dtype('float64')
    z=z.rebin(scale=[4,4,1])
s[0]=s[0].rebin(scale=[4,4,1])
s[1]=s[1].rebin(scale=[4,4,1])
s[2]=s[2].rebin(scale=[4,4,1])    
for z in s:
    z=z.map(gaussian_filter,sigma=3.0)
#%%


s1_avgcounts = s1.inav[:,:].data.sum()/(s1.data.shape[0]*s1.data.shape[1])
s1_totcounts = s1.inav[:,:].data.sum()
#s1=s1.map(gaussian_filter,sigma=3.0)

print('\n')
print("Medelantal counts: "+str(s1_avgcounts))
print("Totalt antal counts: " + str(s1_totcounts))
print('\n')


#%%


#s2 = s2.isig[650.0::]
s2_avgcounts = s2.inav[:,:].data.sum()/(s2.data.shape[0]*s2.data.shape[1])
s2_totcounts = s2.inav[:,:].data.sum()
#s2=s2.map(gaussian_filter,sigma=3.0)
print('\n')
print("Medelantal counts: "+str(s2_avgcounts))
print("Totalt antal counts: " + str(s2_totcounts))
print('\n')


#%%



s3_avgcounts = s3.inav[:,:].data.sum()/(s3.data.shape[0]*s3.data.shape[1])
s3_totcounts = s3.inav[:,:].data.sum()
#s3=s3.map(gaussian_filter,sigma=3.0)
print('\n')
print("Medelantal counts: "+str(s3_avgcounts))
print("Totalt antal counts: " + str(s3_totcounts))
print('\n')

#%%

#hs.plot.plot_spectra([s3.isig[300.0::].inav[::].sum(),s2.isig[300.0::].inav[::].sum()])


sRef.add_elements(['Cu'])
sRef.add_lines(['Cu_Ka'])


im1 = s[0].get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im1,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
plt.title('Line intensites particle 1')


im2 = s[1].get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im2,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
plt.title('Line intensites particle 2')

im3 = s[2].get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im3,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
plt.title('Line intensites particle 3')

#%%
s[0].decomposition()
s[0].plot_explained_variance_ratio()
plt.title('Explained variance particle 1')

s[1].decomposition()
s[1].plot_explained_variance_ratio()
plt.title('Explained variance particle 2')

s[2].decomposition()
s[2].plot_explained_variance_ratio()
plt.title('Explained variance particle 3')


#%%
s[0].decomposition(algorithm = 'NMF', output_dimension = 3)
s1_loads = s[0].get_decomposition_loadings()
s1_facs = s[0].get_decomposition_factors()

s[1].decomposition(algorithm = 'NMF', output_dimension = 3)
s2_loads = s[1].get_decomposition_loadings()
s2_facs = s[1].get_decomposition_factors()

s[2].decomposition(algorithm = 'NMF', output_dimension = 3)
s3_loads = s[2].get_decomposition_loadings()
s3_facs = s[2].get_decomposition_factors()
#%%
hs.plot.plot_images(s1_loads, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            suptitle = 'NMF-loadings particle 1',
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})
plt.title('NMF-loadings particle 1')

hs.plot.plot_images(s2_loads, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            suptitle = 'NMF-loadings particle 2',
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})


hs.plot.plot_images(s3_loads, cmap='mpl_colors',
            axes_decor='off', per_row=2,
            suptitle = 'NMF-loadings particle 3',
            scalebar=[0], scalebar_color='white',
            padding={'top': 0.95, 'bottom': 0.05,
                     'left': 0.05, 'right':0.78})



#%%

for f in s1_facs:
    f.data /= f.data.max()
hs.plot.plot_spectra(s1_facs.isig[:10000.0],style='cascade',padding=-1) 

plt.title('NMF-factors particle 1')
plt.text(x=8631, y=0.6, s='Zn-K$_\\alpha$', color='k')
plt.axvline(8631, c='k', ls=':', lw=0.5)
plt.text(x=1011, y=0.8, s='Zn-L$_\\alpha$', color='k')
plt.axvline(1011, c='k', ls=':', lw=0.5)
plt.axvline(9713, c='k', ls=':', lw=0.5)
plt.text(x=9713, y=0.8, s='Au-L$_\\alpha$', color='k')
plt.axvline(277, c='k', ls=':', lw=0.5)
plt.text(x=2120, y=-0.1, s='Au-M$_\\alpha$', color='k')
plt.axvline(2120, c='k', ls=':', lw=0.5)
plt.text(x=277, y=0.8, s='C-K$_\\alpha$', color='k')

plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
for f in s2_facs:
    f.data /= f.data.max()
hs.plot.plot_spectra(s2_facs.isig[:10000.0],style='cascade',padding=-1) 

plt.title('NMF-factors particle 2')
plt.text(x=8631, y=0.6, s='Zn-K$_\\alpha$', color='k')
plt.axvline(8631, c='k', ls=':', lw=0.5)
plt.text(x=1011, y=0.8, s='Zn-L$_\\alpha$', color='k')
plt.axvline(1011, c='k', ls=':', lw=0.5)
plt.axvline(9713, c='k', ls=':', lw=0.5)
plt.text(x=9713, y=0.8, s='Au-L$_\\alpha$', color='k')
plt.axvline(277, c='k', ls=':', lw=0.5)
plt.text(x=2120, y=-0.1, s='Au-M$_\\alpha$', color='k')
plt.axvline(2120, c='k', ls=':', lw=0.5)
plt.text(x=277, y=0.8, s='C-K$_\\alpha$', color='k')

plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)

for f in s3_facs:
    f.data /= f.data.max()
hs.plot.plot_spectra(s3_facs.isig[:10000.0],style='cascade',padding=-1) 

plt.title('NMF-factors particle 3')
plt.text(x=8631, y=0.6, s='Zn-K$_\\alpha$', color='k')
plt.axvline(8631, c='k', ls=':', lw=0.5)
plt.text(x=1011, y=0.8, s='Zn-L$_\\alpha$', color='k')
plt.axvline(1011, c='k', ls=':', lw=0.5)
plt.axvline(9713, c='k', ls=':', lw=0.5)
plt.text(x=9713, y=0.8, s='Au-L$_\\alpha$', color='k')
plt.axvline(277, c='k', ls=':', lw=0.5)
plt.text(x=277, y=0.8, s='C-K$_\\alpha$', color='k')
plt.text(x=2120, y=-0.1, s='Au-M$_\\alpha$', color='k')
plt.axvline(2120, c='k', ls=':', lw=0.5)
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
#%% 
z = s1_facs.inav[0]

z.set_elements(['Au','Zn'])
z.set_lines(['Au_La','Zn_La'])


kfac = [2.585,1.401]

bw = z.estimate_background_windows(line_width=[1.7, 2.0])
z.plot(background_windows=bw)

intensities = z.get_lines_intensity(background_windows=bw)#
atomic_percent = z.quantification(intensities, method='CL',composition_units='weight',
                                  factors=kfac)


print('Au at%: ' + str(atomic_percent[0].data[0]))
print('Zn at% ' + str(atomic_percent[1].data[0]))


