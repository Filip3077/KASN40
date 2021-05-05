# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 14:20:19 2021

@author: Lucas
"""
import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from coreshellfunctions import*
#%%

sAg = hs.load("../Spectra/20nm cube Cu0Ag100.msa",signal_type="EDS_TEM")
sCu = hs.load("../Spectra/20nm cube Cu100Ag0.msa",signal_type="EDS_TEM")
sC=hs.load("../Spectra/Carbonbackground.msa", signal_type="EDS_TEM")
cal = sAg #0.45*sAg.data+0.45*sCu.data+0.1*sC.data
specCore = 0.7*sCu.data + 0.3*sAg.data; specShell = 0.9*sAg.data + 0.1*sCu.data

intensity = 1/20
a = CoreShellP(50, 20, 15, intensity, intensity, 1)
a = CoreShellSpec(a, specCore, specShell)
a.add_background(sC,1)
p = a.core+a.shell; core = a.core; shell = a.shell
p.add_poissonian_noise()
p = setCalibration(p, cal); core = setCalibration(core, cal); shell = setCalibration(shell, cal)
cut_spectrum_bottom(p, 500); cut_spectrum_bottom(core, 500); cut_spectrum_bottom(shell, 500)

sumsignal = np.sum(p.data, axis=(0,1,2))
meansignal = np.sum(p.data, axis=(0,1,2))/50**2

print("Total counts: "+str(sumsignal))
print("Mean counts: "+str(meansignal))

newbin = [1, 1, 8]
p.change_dtype('float64'); core.change_dtype('float64'); shell.change_dtype('float64')
p = p.rebin(scale = newbin); core = core.rebin(scale = newbin); shell = shell.rebin(scale = newbin)

#%%

im = p.get_lines_intensity()

hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
    colorbar='single', vmin='1th', vmax='99th', scalebar='all',
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
             'right':0.85, 'wspace':0.20, 'hspace':0.10})


#%%

nfac = 2 #int(input("Chose nr of components :"))
p.decomposition(True, algorithm='NMF', output_dimension=nfac)

loadings = p.get_decomposition_loadings().inav[0:nfac]
factors = p.get_decomposition_factors().inav[0:nfac]

hs.plot.plot_spectra(factors,style='cascade')
plt.title('NMF')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    
hs.plot.plot_images(loadings, cmap='mpl_colors', suptitle='NMF',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%

def RotCompactness(a, *loads):

    return np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1]).sum(axis=(0,1))

def RotCompactnessMoment(a, *loads):
    #Gör en rotation och ger "vridmomentet" för kärnarn. För en "kompakt"
    #kärna ska detta vara litet eftersom all intensitet är nära mitten.
    size = len(loads[-1])
    ld = np.abs(np.linspace(0,size-1,size).reshape((size,1))@np.ones(size).reshape((1,size))-(size-1)/2)
    ld = np.sqrt(ld**2 + ld.T**2)
    return np.sum(np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1])*ld**4, axis=(0,1))

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

def CorrectShell(b, *args):
    temp_shell = args[1]+b*args[0]
    temp_profile = radial_profile(temp_shell, [24.5, 24.5])
    
    return np.sum((temp_profile/np.max(temp_profile)-args[2])**2)

def transfer_elements(factors,loadings,c,s):
    
    size = loadings.inav[-1].data.shape[0]

    res = minimize_scalar(RotCompactnessMoment, args=tuple(loadings.data), bounds=(-np.pi/2, np.pi/2), method='bounded')

    a = res.x
    loadings_rot = loadings.deepcopy(); factors_rot = factors.deepcopy()
    loadings_rot.data = [np.cos(a)*loadings.inav[c].data+np.sin(a)*loadings.inav[s].data, (np.cos(a)*loadings.inav[s].data-np.sin(a)*loadings.inav[c].data) ]
    factors_rot.data = [np.cos(a)*factors.inav[c].data+np.sin(a)*factors.inav[s].data, (np.cos(a)*factors.inav[s].data-np.sin(a)*factors.inav[c].data) ]
    
    dest_tot = np.sqrt(np.sum((loadings_rot.inav[c].data+loadings_rot.inav[s].data >0.01), axis=(0,1))/np.pi)
    dest_core = np.sqrt(np.sum((loadings_rot.inav[c].data >0.01), axis=(0,1))/np.pi)
    dest_shell = dest_tot-dest_core
    print('dest_shell is' + str(dest_shell))

    # center_int_frac = dest_shell/np.sqrt(dest_tot**2-dest_core**2)

    shell_profile = radial_profile(loadings_rot.inav[s].data, [(size-1)/2, (size-1)/2])
    shell_profile_ideal = np.arange(0,shell_profile.shape[0],1)**2
    shell_profile_ideal = np.nan_to_num(np.sqrt(dest_tot**2-shell_profile_ideal)) - np.nan_to_num(np.sqrt(dest_core**2-shell_profile_ideal))
    shell_profile_ideal= shell_profile_ideal / np.max(shell_profile_ideal)
    loadings_tup=[loadings_rot.inav[c].data,loadings_rot.inav[s].data]
    res1 = minimize_scalar(CorrectShell, args=tuple(loadings_tup)+(shell_profile_ideal,), bounds=(-0.8, 0.8), method='bounded')
    b = res1.x
    print('B is '+str(b))
    loadings_corr=loadings_rot.deepcopy(); factors_corr=factors_rot.deepcopy()
    loadings_corr.inav[s].data=loadings_rot.inav[s].data+b*loadings_rot.inav[c].data
    # loadings_corr.inav[c].data=loadings_rot.inav[c].data-b*loadings_rot.inav[s].data
    factors_corr.inav[c] = factors_rot.inav[c].data-b*factors_rot.inav[s].data
    return factors_corr,loadings_corr

factors_corr, loadings_corr = transfer_elements(factors, loadings, 0,1)



#%%

hs.plot.plot_spectra(factors_corr,style='cascade')
plt.title('NMF_Corr')
plt.text(x=8040, y=0.8, s='Cu-K$_\\alpha$', color='k')
plt.axvline(8040, c='k', ls=':', lw=0.5)
plt.text(x=930, y=0.8, s='Cu-L$_\\alpha$', color='k')
plt.axvline(930, c='k', ls=':', lw=0.5)
plt.axvline(2984, c='k', ls=':', lw=0.5)
plt.text(x=2984, y=0.8, s='Ag-L$_\\alpha$', color='k')

    
hs.plot.plot_images(loadings_corr, cmap='mpl_colors', suptitle='NMF_Corr',
                    axes_decor='off', per_row=3,
                    scalebar=[0], scalebar_color='white',
                    padding={'top': 0.95, 'bottom': 0.05,
                             'left': 0.05, 'right':0.78})

#%%

def quantify(spec, kfac=[1,0.72980399], linewidth=[5.0,7.0], stack=False):
    if stack:
        for s in spec:
            print(s)
            bw = s.estimate_background_windows(line_width=linewidth)
            intensities = s.get_lines_intensity(background_windows=bw)
            wtAg = s.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data[0]
            wtCu = s.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data[0]
    else:
        bw = spec.estimate_background_windows(line_width=linewidth)
        intensities = spec.get_lines_intensity(background_windows=bw)
        wtAg = spec.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[0].data[0]
        wtCu = spec.quantification(intensities, method='CL', factors=kfac,composition_units='weight')[1].data[0]
    
    return wtAg, wtCu

NMF_p = cLoadsFacs(loadings, factors).split()
NMF_corr_p =cLoadsFacs(loadings_corr, factors_corr).split()
# NMF_corr_p[0].data = NMF_corr_p[0].data-np.min(NMF_corr_p[0].data)
# NMF_corr_p[1].data = NMF_corr_p[1].data-np.min(NMF_corr_p[1].data)

core_im = core.get_lines_intensity(); shell_im = shell.get_lines_intensity()
core_shell_im = [core_im[0],core_im[1],shell_im[0],shell_im[1]]

temp = hs.signals.BaseSignal(np.random.random((50,50)))
NMF_im = [temp,temp,temp,temp]
NMF_corr_im = [temp,temp,temp,temp]
wtp = []; wtp.append(quantify(core.inav[:,:].sum())); wtp.append(quantify(shell.inav[:,:].sum()))

z = 0
for i in range(2):
    NMF_p[i] = setCalibration(NMF_p[i], p.inav[24,24] )
    NMF_corr_p[i] = setCalibration(NMF_corr_p[i], p.inav[24,24])
    wtp.append(quantify(NMF_p[i].inav[:,:].sum())); wtp.append(quantify(NMF_corr_p[i].inav[:,:].sum()))
    pc = SpecErrAbs2D(NMF_p[i], core)
    ps = SpecErrAbs2D(NMF_p[i], shell)
    cc = SpecErrAbs2D(NMF_corr_p[i], core)
    cs = SpecErrAbs2D(NMF_corr_p[i], shell)
    print('NMF PC: '+str(i)+' core match: '+str(pc))
    print('NMF PC: '+str(i)+' shell match: '+str(ps))
    print('NMF_Corr PC: '+str(i)+' core match: '+str(cc))
    print('NMF_Corr PC: '+str(i)+' shell match: '+str(cs))
    temp_im = NMF_p[i].get_lines_intensity()
    temp_corr_im = NMF_corr_p[i].get_lines_intensity()
    for j in range(2):
        NMF_im[j+z] = temp_im[j]
        NMF_corr_im[j+z] = temp_corr_im[j]
    z = 2


hs.plot.plot_images(core_shell_im, cmap='RdYlBu_r', axes_decor='off', #tight_layout=True,
    colorbar='single', vmin='1th', vmax='99th', scalebar='all', per_row=2,
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

hs.plot.plot_images(NMF_im, cmap='RdYlBu_r', axes_decor='off', #tight_layout=True,
    colorbar='single', vmin='1th', vmax='99th', scalebar='all', per_row=2,
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

hs.plot.plot_images(NMF_corr_im, cmap='RdYlBu_r', axes_decor='off', #tight_layout=True
    colorbar='single', vmin='1th', vmax='99th', scalebar='all', per_row=2,
    scalebar_color='black', suptitle_fontsize=16,
    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
              'right':0.85, 'wspace':0.20, 'hspace':0.10})

#%%

kfac = [1,0.72980399]

bw = p.inav[:,:].sum().estimate_background_windows(line_width=[5.0, 7.0])

p.inav[:,:].sum().plot(background_windows=bw)

p.inav[:,:].sum().get_lines_intensity(background_windows=bw, plot_result=True)






#%%
# loadings = loadings.data; factors = factors.data
# from scipy.optimize import minimize_scalar

# def checkshell(loadings):
#     #Checkshell ser bara till att 2:a PC har positiv intensitet i skalet
#     x, y = loadings[0].shape
#     xc = int(x/2)
#     yc = int(y/2)
    
#     return -np.mean(loadings[1][yc-2:yc+2, xc-2:xc+2])/np.abs(np.mean(loadings[1][yc-2:yc+2, xc-2:xc+2]))

# def RotCompactnessMoment(a, *loads):
#     #Gör en rotation och ger "vridmomentet" för kärnarn. För en "kompakt"
#     #kärna ska detta vara litet eftersom all intensitet är nära mitten.
#     return np.sum(np.abs(np.cos(a)*loads[0]+np.sin(a)*loads[1])*ld**4, axis=(0,1))

# def radial_profile(data, center):
#     y, x = np.indices((data.shape))
#     r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
#     r = r.astype(np.int)

#     tbin = np.bincount(r.ravel(), data.ravel())
#     nr = np.bincount(r.ravel())
#     radialprofile = tbin / nr
#     return radialprofile 

# def CorrectShell(b, *args):
#     temp_shell = args[1]+b*args[0]
#     temp_profile = radial_profile(temp_shell, [24.5, 24.5])
    
#     return np.sum((temp_profile/np.max(temp_profile)-args[2])**2)

# #Test så att principalkomponent 2 har positiv loading i skalet
# t = checkshell(loadings)
# loadings[1] = loadings[1]*t
# factors[1] = factors[1]*t

# keV = np.linspace(-0.2, -0.2+0.01*2048, int(2048/newbin[2]))

# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
# ax1.imshow(loadings[1])
# ax2.plot(keV, factors[1])
# ax2.set_xlim([-0.1, 10])

# #Rotering för att få fram en ren kärn-loading
# res = minimize_scalar(RotCompactnessMoment, args=tuple(loadings), bounds=(-np.pi/2, np.pi/2), method='bounded')

# a = res.x

# loadings_rot = [np.cos(a)*loadings[0]+np.sin(a)*loadings[1], (np.cos(a)*loadings[1]-np.sin(a)*loadings[0]) ]
# factors_rot = [np.cos(a)*factors[0]+np.sin(a)*factors[1], (np.cos(a)*factors[1]-np.sin(a)*factors[0]) ]

# #Uppskatta storlek på kärna och kärna+skal utifrån deras areor. Enklast möjliga
# #metod. Sist skapas en "ideal" radiell profil som visar hur skalets loading
# #bör se ut.
# dest_tot = np.sqrt(np.sum((loadings[0] >0.01), axis=(0,1))/np.pi)
# dest_core = np.sqrt(np.sum((loadings_rot[0] >0.01), axis=(0,1))/np.pi)
# dest_shell = dest_tot-dest_core

# center_int_frac = dest_shell/np.sqrt(dest_tot**2-dest_core**2)

# shell_profile = radial_profile(loadings_rot[1], [24.5, 24.5])

# shell_profile_ideal = np.arange(0,shell_profile.shape[0],1)**2
# shell_profile_ideal = np.nan_to_num(np.sqrt(dest_tot**2-shell_profile_ideal)) - np.nan_to_num(np.sqrt(dest_core**2-shell_profile_ideal))
# shell_profile_ideal= shell_profile_ideal / np.max(shell_profile_ideal)

# #Skalet har vid det här laget oftast antingen ett "hål" i mitten, eller är 
# #tjockare uppe på partikeln. Kärnans loading läggs till eller dras ifrån skalet
# #för att skapa en skal-loading som har "rätt" profil.
# res1 = minimize_scalar(CorrectShell, args=tuple(loadings_rot)+(shell_profile_ideal,), bounds=(-0.8, 0.8), method='bounded')

# b = res1.x

# #Skalets loading korrigeras för att ge rätt profil. Ändringen kompenseras 
# #genom att ändra kärnans factor med motsvarande konstant: ex. om skalet hade 
# #"hål" från början måste lite silver i detta fall läggas till kärnan för att 
# #skalet ska kunna täcka kärnan jämnt.
# loadings_corr = [loadings_rot[0], loadings_rot[1]+b*loadings_rot[0]]
# factors_corr = [factors_rot[0]-b*factors_rot[1], factors_rot[1] ]

# shell_profile2 = radial_profile(loadings_corr[1], [24.5, 24.5])
# core_profile = radial_profile(loadings_corr[0], [24.5, 24.5])

#%%

# def radial_profile(data, center):
#     y, x = np.indices((data.shape))
#     r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
#     r = r.astype(np.int)

#     tbin = np.bincount(r.ravel(), data.ravel())
#     nr = np.bincount(r.ravel())
#     radialprofile = tbin / nr
#     return radialprofile 

# def CorrectShell(b, *args):
#     temp_shell = args[1]+b*args[0]
#     temp_profile = radial_profile(temp_shell, [24.5, 24.5])
    
#     return np.sum((temp_profile/np.max(temp_profile)-args[2])**2)

# def transfer_elements(factors,loadings,c,s,size):
#     dest_tot = np.sqrt(np.sum((loadings.inav[c].data+loadings.inav[s].data >0.01), axis=(0,1))/np.pi)
#     print('Dest_tot'+str(dest_tot))
#     dest_core = np.sqrt(np.sum((loadings.inav[c].data >0.01), axis=(0,1))/np.pi)
#     print('Dest_core'+str(dest_core))
#     # dest_core = 15
    
#     dest_shell = dest_tot-dest_core
#     print('Dest_shell'+str(dest_shell))


#     center_int_frac = dest_shell/np.sqrt(dest_tot**2-dest_core**2)

#     shell_profile = radial_profile(loadings.inav[s].data, [(size-1)/2, (size-1)/2])
#     shell_profile_ideal = np.arange(0,shell_profile.shape[0],1)**2
#     print('Shell_profile_ideal is '+str(shell_profile_ideal))
#     shell_profile_ideal = np.nan_to_num(np.sqrt(dest_tot**2-shell_profile_ideal)) - np.nan_to_num(np.sqrt(dest_core**2-shell_profile_ideal))
#     shell_profile_ideal= shell_profile_ideal / np.max(shell_profile_ideal)
#     loadings_tup=[loadings.inav[c].data,loadings.inav[s].data]
#     res1 = minimize_scalar(CorrectShell, args=tuple(loadings_tup)+(shell_profile_ideal,), bounds=(-0.8, 0.8), method='bounded')
#     b = res1.x
#     print('B is '+str(b))
#     loadings_corr=loadings.deepcopy()
#     factors_corr=factors.deepcopy()
#     fcorr = [factors.inav[c].data-b*factors.inav[s].data, factors.inav[s].data+b*factors.inav[c].data]
#     lcorr = [loadings.inav[c].data-b*loadings.inav[s].data, loadings.inav[s].data+b*loadings.inav[c].data]
#     # loadings_corr.inav[s].data=loadings.inav[s].data+(-0.4)*loadings.inav[c].data
#     # factors_corr.inav[c].data=factors.inav[c].data-(-0.4)*factors.inav[s].data
#     factors_corr.data = fcorr; loadings_corr.data = lcorr
#     return factors_corr,loadings_corr

