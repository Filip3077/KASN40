# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 08:17:19 2021

@author: Filip
"""

import hyperspy.api as hs
import matplotlib.pyplot as plt


#s20nm =  hs.load("EDS Data_particle*.rlp",stack=True,signal_type="EDS_TEM")

s1 = hs.load("AuZn4.rpl",stack=True,signal_type="EDS_TEM").T
s2 = hs.load("AuZn2.rpl",stack=True,signal_type="EDS_TEM").T
s3 = hs.load("AuZn3.rpl",stack=True,signal_type="EDS_TEM").T
sRef =  hs.load("PureCuFilip.msa",stack=True,signal_type="EDS_TEM")
sRef = sRef.isig[0:1500]


s1.axes_manager[0].name = 'y'
s1.axes_manager[1].name = 'x'
s1.axes_manager[-1].name = 'E'
s1.axes_manager.signal_axes[0].units = 'eV'
s1.axes_manager['E'].scale = 10
s1.axes_manager['x'].units = 'nm'
s1.axes_manager['y'].units = 'nm'
s1.axes_manager['E'].offset = -200
s1.axes_manager['y'].scale = 0.22
s1.axes_manager['x'].scale = 0.22
#s1.isig[200.0::].inav[::].sum().plot()
s1 = s1.rebin(scale = [4,4,1])

s1_avgcounts = s1.inav[:,:].data.sum()/(s1.data.shape[0]*s1.data.shape[1])
s1_totcounts = s1.inav[:,:].data.sum()

print('\n')
print("Medelantal counts: "+str(s1_avgcounts))
print("Totalt antal counts: " + str(s1_totcounts))
print('\n')




s2.axes_manager[0].name = 'y'
s2.axes_manager[1].name = 'x'
s2.axes_manager[-1].name = 'E'
s2.axes_manager.signal_axes[0].units = 'eV'
s2.axes_manager['E'].scale = 10
s2.axes_manager['x'].units = 'nm'
s2.axes_manager['y'].units = 'nm'
s2.axes_manager['E'].offset = -200
s2.axes_manager['y'].scale = 0.22
s2.axes_manager['x'].scale = 0.22
#s2.isig[200.0::].inav[::].sum().plot()
s2 = s2.rebin(scale = [4,4,1])
#s2 = s2.isig[650.0::]
s2_avgcounts = s2.inav[:,:].data.sum()/(s2.data.shape[0]*s2.data.shape[1])
s2_totcounts = s2.inav[:,:].data.sum()

print('\n')
print("Medelantal counts: "+str(s2_avgcounts))
print("Totalt antal counts: " + str(s2_totcounts))
print('\n')




s3.axes_manager[0].name = 'y'
s3.axes_manager[1].name = 'x'
s3.axes_manager[-1].name = 'E'
s3.axes_manager['x'].units = 'nm'
s3.axes_manager['y'].units = 'nm'
s3.axes_manager.signal_axes[0].units = 'eV'
s3.axes_manager['E'].scale = 10
s3.axes_manager['E'].offset = -200
s3.axes_manager['y'].scale = 0.22
s3.axes_manager['x'].scale = 0.22
s3.isig[200.0::].inav[::].sum().plot()
s3 = s3.rebin(scale = [4,4,1])

s3_avgcounts = s3.inav[:,:].data.sum()/(s3.data.shape[0]*s3.data.shape[1])
s3_totcounts = s3.inav[:,:].data.sum()

print('\n')
print("Medelantal counts: "+str(s3_avgcounts))
print("Totalt antal counts: " + str(s3_totcounts))
print('\n')


#hs.plot.plot_spectra([s3.isig[300.0::].inav[::].sum(),s2.isig[300.0::].inav[::].sum()])
s1.add_elements(['Au','Zn','C','C']) #Lägger in element igen tydligen förs de inte med 
s1.add_lines(['Au_La','Zn_La','C_Ka','C_Ka'])

s2.add_elements(['Au','Zn','C']) #Lägger in element igen tydligen förs de inte med 
s2.add_lines(['Au_La','Zn_La','C_Ka'])


s3.add_elements(['Au','Zn','C']) #Lägger in element igen tydligen förs de inte med 
s3.add_lines(['Au_La','Zn_La','C_Ka'])

sRef.add_elements(['Cu'])
sRef.add_lines(['Cu_Ka'])


im1 = s1.get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im1,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
plt.title('Line intensites particle 1')


im2 = s2.get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im2,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
plt.title('Line intensites particle 2')

im3 = s3.get_lines_intensity()  #Tar ut två "bilder" en för varje xray signal dvs Ag och Cu
hs.plot.plot_images(im3,  cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})
plt.title('Line intensites particle 3')

#%%
s1.decomposition()
s1.plot_explained_variance_ratio()
plt.title('Explained variance particle 1')

s2.decomposition()
s2.plot_explained_variance_ratio()
plt.title('Explained variance particle 2')

s3.decomposition()
s3.plot_explained_variance_ratio()
plt.title('Explained variance particle 3')


#%%
s1.decomposition(algorithm = 'NMF', output_dimension = 3)
s1_loads = s1.get_decomposition_loadings()
s1_facs = s1.get_decomposition_factors()

s2.decomposition(algorithm = 'NMF', output_dimension = 3)
s2_loads = s2.get_decomposition_loadings()
s2_facs = s2.get_decomposition_factors()

s3.decomposition(algorithm = 'NMF', output_dimension = 3)
s3_loads = s3.get_decomposition_loadings()
s3_facs = s3.get_decomposition_factors()
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
s = s3.sum().sum()#s3_facs.inav[0]

s.set_elements(['Au','Zn'])
s.set_lines(['Au_La','Zn_La'])


kfac = [2.585,1.401]

bw = s.estimate_background_windows(line_width=[1.7, 2.0])
s.plot(background_windows=bw)

intensities = s.get_lines_intensity(background_windows=bw)
atomic_percent = s.quantification(intensities, method='CL',
                                  factors=kfac)


print('Au at%: ' + str(atomic_percent[0].data[0]))
print('Zn at% ' + str(atomic_percent[1].data[0]))


