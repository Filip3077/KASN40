# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:55:56 2021

@author: Lucas
"""

import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from numpy import linalg
import math
import scipy.misc


#%% Proc/loadassign.py
"""
Created on Fri Mar  5 14:38:15 2021

@author: Jonas
"""



def setCalibration(ucMap,calSpec,background=False):
    ucMap.set_signal_type("EDS_TEM")
    ucMap.axes_manager[0].name = 'y'
    ucMap.axes_manager[1].name = 'x'
    ucMap.axes_manager['x'].units = 'nm'
    ucMap.axes_manager['y'].units = 'nm'
    ucMap.axes_manager[-1].name = 'E'
    ucMap.get_calibration_from(calSpec)
    if background:
        ucMap.add_elements(['Ag','Cu', 'C'])
        ucMap.add_lines(['Ag_La','Cu_Ka','C_Ka'])
    else:
        ucMap.add_elements(['Ag','Cu'])
        ucMap.add_lines(['Ag_La','Cu_Ka'])
    return ucMap

def AuZnCalibration(im,pixelScale=1/6,elements=['Au','Zn'],peaks=['Au_La','Zn_La']):
    im.set_signal_type("EDS_TEM")
    im.axes_manager[0].name = 'y'
    im.axes_manager[1].name = 'x'
    im.axes_manager['x'].units = 'nm'
    im.axes_manager['y'].units = 'nm'
    im.axes_manager['x'].scale = pixelScale
    im.axes_manager['y'].scale = pixelScale
    im.axes_manager[-1].name = 'E'
    im.axes_manager['E'].units = 'eV'
    im.axes_manager['E'].scale = 10
    im.axes_manager['E'].offset=-200
    im.add_elements(elements)
    im.add_lines(peaks)
    return im

def redbluePlot(im,Title=''):
    im = im.get_lines_intensity()

    hs.plot.plot_images(im, suptitle=Title, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
        colorbar='single', vmin='1th', vmax='99th', scalebar='all',
        scalebar_color='black', suptitle_fontsize=16,
        padding={'top':0.8, 'bottom':0.10, 'left':0.05,
                 'right':0.85, 'wspace':0.20, 'hspace':0.10}) 
    return 0



def cLoadsFacs(loads,facs):
    #Antar att både loads och facs kommer från samma "ursprung" och har samma ordning och dimentioner
    #För att få ett korrekta dimentioner på  hyperspy objektet böhöver loads transponeras från [| x y]  till [x y |] har att göra med hur energiaxeln behandlas 
    
    dim = len(loads)
    size = len(loads.isig)
    esize = len(facs.isig)
    combinedMat = np.empty((dim,size,size,esize))
    
    for i in range(dim):
        #För att få ett korrekta dimentioner på  hyperspy objektet behöver loads transponeras från [| x y]  till [x y |] har att göra med hur energiaxeln behandlas 
        combinedMat[i] = (loads.inav[i].T*facs.inav[i]).data
    
    
    combined = hs.signals.BaseSignal(combinedMat)
    combined=combined.transpose(signal_axes=[0],navigation_axes=[1,2,3])
    s = loads.metadata.General.title.split("of")
    combined.metadata.General.title = "Reconstructed: " + s[1]
    return combined

def checkLoadOnly(core,shell,statload,components=2,method='abs'):
   statcoretest = None
   statshelltest = None
   statcore = 1000
   statshell = 1000
   if method=='abs':
       for i in range(components):
           statcoretest=SpecErrAbs2D(statload.inav[i].data, core,signal=False)
           statshelltest=SpecErrAbs2D(statload.inav[i].data, shell,signal=False)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
                statshell = statshelltest
                index_s = i
           if index_c == index_s:
               print('Warning! It guessed the same component for core and shell.')
   elif method=='neuc':
       print("neuc hasn't been implemented for this function yet!")
       for i in range(components):
           statcoretest=SpecErrNEuc(statload.inav[i], core)
           statshelltest=SpecErrNEuc(statload.inav[i], shell)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
               statshell = statshelltest
               index_s = i
           if index_c == index_s:
              print('Warning! It guessed the same component for core and shell.')
   else:
       print(method+" is not a valid comparison method. Enter 'abs' or 'neuc'")
   return index_c,index_s #index för den som bäst passar core respektive shell

def checkLoadFit(core,shell,statfac,statload,components=2,method='abs'):
    
   '''Takes the core and shell of a core-shell particle (Hyperspy-signals) as well as \n
   the factors and loadings of a statistical method and returns the indices\n
   of the factor+loading combo that best correspond to core and shell respectively.
    '''
   Statspec = cLoadsFacs(statload, statfac)
   statcoretest = None
   statshelltest = None
   statcore = np.inf
   statshell = np.inf
   if method=='abs':
       for i in range(components):
           statcoretest=SpecErrAbs2D(Statspec.inav[i], core)
           statshelltest=SpecErrAbs2D(Statspec.inav[i], shell)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
                statshell = statshelltest
                index_s = i
           if index_c == index_s:
               print('Warning! It guessed the same component for core and shell.')
   elif method=='neuc':
        for i in range(components):
           statcoretest=SpecErrNEuc(Statspec.inav[i], core)
           statshelltest=SpecErrNEuc(Statspec.inav[i], shell)
           if statcoretest < statcore:
               statcore = statcoretest
               index_c = i
           if statshelltest < statshell:
               statshell = statshelltest
               index_s = i
        if index_c == index_s:
              print('Warning! It guessed the same component for core and shell.')
   else:
       print(method+" is not a valid comparison method. Enter 'abs' or 'neuc'")
   return index_c,index_s #index för den som bäst passar core respektive shell

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

#%% Proc/specerr.py
"""
Created on Mon Feb  8 13:18:03 2021

@author: Jonas
"""

def SpecErrAbs (spec1, spec2):
    '''Compares two spectra as vectors and generates a fractional difference between them as \n
    a percentage error of spec1'''
    l1=len(spec1);
    l2=len(spec2);
    s=0;
    if (l1>=l2):# If the spectra have different lengths the longest is cut to match the shorter one.  
        l=l2;
    else:
        l=l1;
    for i in range(l):
        s=s+abs(spec1[i]-spec2[i]);
    if (sum(spec1)==0):
        return 0
    else:
        return s/sum(spec1);

def SpecErrSq (spec1, spec2):
    ''' Takes two spectra as vectors and returns the coefficient of determination using spec2 as a model of spec1'''
    l1=len(spec1);
    l2=len(spec2);
    s=0;
    z=0;
    if (l1>=l2):
           l=l2;
    else:
           l=l1;
    for i in range(l):
            s=s+((spec1[i]-spec2[i]))**2;
            z=z+((spec1[i]-np.mean(spec1)))**2            
    return 1-s/z;# Returns R^2 where spec2 is seen as a model of spec1.

def SpecErrAbs2D(map1,relMap,signal=True):
    '''Takes a numpy array with spectra stored in each element and generates\n
    a fractional difference between the two arrays based in the spectra. The norm \n
    is spec1.'''
    if signal:
        reldif = abs(map1-relMap).data.sum()/relMap.data.sum() 
        return reldif
    else:
        l1=relMap.shape;
        l2=map1.shape;
        k0=min(l1[0],l2[0]);
        k1=min(l1[1],l2[1]);
        store=np.empty((k0,k1))
        for i in range(k0):
            for j in range(k1):
                store[i][j]=abs(relMap[i][j]-map1[i][j])
        return np.sum(store)/np.sum(relMap)

def SpecErrNEuc(map1, relMap):
    diffmap=(map1-relMap)**2;
    reldif=diffmap.data.sum()/relMap.data.sum()
    return reldif;

#%% Proc/specMapDiff.py
"""
Created on Fri Feb 19 15:39:07 2021

@author: Filip
"""

def specMapDiff(map1,refMap):
    #Om map1 och map2 är hyperspy objekt går det helt enkelt att ta differensen direkt samt att ta absolutvärdet av denna. Om dimentionerna stämmer dvs. 
    diff = abs(map1-refMap)
    return diff


def rel(EF_map,ref):
    '''
    Tar den relativa 
    '''
    refMap = EF_map.data/ref.data
    where_are_NaNs = np.isnan(refMap)
    refMap[where_are_NaNs] = 0
    refMap = hs.signals.BaseSignal(refMap).T
    return refMap

#%% rot/radial_profile.py
"""
Created on Wed Apr  7 16:49:39 2021

@author: Jonas
"""

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
    # loadings_corr.inav[s].data=loadings_rot.inav[s].data+b*loadings_rot.inav[c].data
    # loadings_corr.inav[c].data=loadings_rot.inav[c].data-b*loadings_rot.inav[s].data
    # factors_corr.inav[c] = factors_rot.inav[c].data-b*factors_rot.inav[s].data
    return factors_corr,loadings_corr

#%% rot/varimax.py
"""
Created on Mon Nov  2 21:39:32 2020

@author: Martin
"""

def varimax(Phi , gamma = 1.0 , q = 100  , tol = 1e-10):
    ''' The original Varimax function given by Martin\n
    Takes the unfolded (2+1D -> 1+1D) loadings matrix as main variable. '''
    import numpy as np
    from numpy import linalg
    p, k = Phi.shape
    R = np.eye(k)
    d=1e-6
    
    for i in range(q) :
        d_old = d
        Lambda = np.dot( Phi , R)
        temp1 = np.dot( Lambda , np.diag( np.diag( np.dot( Lambda.T , Lambda ) ) ) )
        temp2 = np.dot( Phi.T , np.asarray( Lambda )**3 - ( gamma/ p ) * temp1 )
        u,s,vh = linalg.svd(temp2)
        R = np.dot ( u , vh )
        d = np.sum( s )
        
        # print(d/d_old)
        
        if d/d_old < (1+tol) :
            break
    return R

def FullVarimax(factors,loadings,components=-1):
    ''' Takes factors and loadings and returns (positive) rotated factors and loadings\n
    Optional variable components should be used if the number of components have not been choosen before.'''
    import numpy as np
    from numpy import linalg
    
    if components == -1:
        nfac = len(factors)
    else:
        nfac = components
        
    factors_selected = factors.inav[0:nfac]
    loadings_selected = loadings.inav[0:nfac]

    #Unfold to turn loadings 2+1D matrix into 1+1D matrix for varimaxfunction.
    loadings_unfold = loadings_selected.deepcopy()
    loadings_unfold.unfold()
    factors_temp = factors_selected.deepcopy()
    factors_selected = factors_selected.data
    loadings_selected = loadings_unfold.data

    R =varimax(loadings_selected.T)

    loadings_selected_rot = np.matmul(loadings_selected.T, R).T
    factors_selected_rot = np.matmul(linalg.inv(R), factors_selected)
    
    #Flipping negative factors and loadings. Only guaranteed for fully negative (non-mixed negative and positive peaks)
    for i in range(nfac):
        if factors_selected_rot[i,:].sum() < 0:
            factors_selected_rot[i,:] = -factors_selected_rot[i,:]
            if loadings_selected_rot[i,:].sum() < 0:
                loadings_selected_rot[i,:] = -loadings_selected_rot[i,:]
    
    factors_temp.data = factors_selected_rot # this way factors is kept as hs objekt.
    #Fold the loadings matrix to original format.
    loadings_unfold.data = loadings_selected_rot
    loadings_unfold.fold()
    
    return factors_temp, loadings_unfold

#%% sim/coreshellp.py
"""
Created on Tue Feb  9 14:20:11 2021

@author: Jonas
"""

def cut_spectrum_bottom(hz,at):
    at=float(at);
    hz.isig[0:at].data-=hz.isig[0:at].data;
    return None
def cut_spectrum_top(hz,at):
    at=float(at);
    hz.isig[at:-1].data-=hz.isig[at:-1].data;
    return None
def cut_spectrum_range(hz,start,stop):
    start,stop=float(start),float(stop);
    hz.isig[start:stop].data-=hz.isig[start:stop];
    return None

class EdxMat:
    def __init__(self,size,r: float,dens:float,l: float):
        '''size: defines the 'picture' size as size x size pixels\n
        r: defines the radius of the spherical particle\n
        dens: density of the material the particle is made of\n
        l: pixel size in the same unit as r (i.e. the area is l^2)'''
        self.size=size;
        self.r=r;
        self.dens=dens;
        self.l=l;
        rred=r/l;
        self.mid= float((size-1)/2);#-1 because indices start from zero
        #The constructor now constructs a matrix using the input arguments
        mat=np.zeros((size,size));
        self.mat=mat;
        for i in range(size):
            for j in range(size):
                if rred**2>=((i-self.mid)**2+(j-self.mid)**2):
                    mat[i][j]=2*dens*self.__thick(i,j);
        
    def __thick(self,n,m):
        nr=n-self.mid;
        mr=m-self.mid;
        l=self.l;
        r=self.r;
        return math.sqrt(r**2-(nr)**2-(mr)**2)*l**2; 
    

class CoreShellP:
    def __init__(self,size,r1:float,r2:float,dens1:float,dens2:float,l:float):
        '''Constructs an object modeling a core-shell spherical particle.\n
        The object properties core and shell  are matrices \n
        modeling the core and shell respectively.\n
        All arguments are defined as for EdxMat.'''
        
        x=EdxMat(size,r1,dens1,l);
        y=EdxMat(size,r2,dens2,l);
        self.size=size;
        self.l=l;
        if r1>=r2:
            z=EdxMat(size,r2,dens1,l)
            self.shell=x.mat-z.mat;
            self.core=y.mat;

        else:
            z=EdxMat(size,r1,dens2,l)
            self.shell=y.mat-z.mat;
            self.core=x.mat;
            


class CoreShellSpec:
     def __init__(self,a,specCore,specShell,signal=True):
         '''Takes a CoreShellP object a and inserts the core spectrum \n
         spec1 and the shell spectrum spec2 at the appropriate places.\n
         Both spectra must be HyperSpy 1D signals.'''
         spec1 = specCore; spec2 = specShell
         self.signal=signal
         L=len(spec1.data)
         arr=np.empty((a.size,a.size,L))
         core=arr.copy()
         shell = arr.copy()
         for i in range(0,a.size):
             for j in range(0,a.size):
                 core[i,j,0:L]=a.core[i,j]*spec1.data
                 shell[i,j,0:L]=a.shell[i,j]*spec2.data
         if signal:
             self.core = hs.signals.Signal1D(core)
             self.shell = hs.signals.Signal1D(shell)
         else:
             self.core=core;
             self.shell=shell;
         self.full=self.core+self.shell   #Ensures backwards compatibility
         self.base=a;
     def getmatr(self):
         '''Returns the matrix representing the entire core-shell particle'''
         return self.core+self.shell
     def is_signal(self):
         return self.signal;
     def add_background(self,backspec, thickness):
         '''Adds a background to the core-shell particle. The thickness is a multiplier\n
         on the original intensity.'''
         for i in range(self.base.size):
             for j in range(self.base.size):
                 if self.signal:
                     self.core.data[i][j]+=backspec.data*thickness;
                 else:
                     self.core[i][j]+=backspec.data*thickness;
                     
     def add_core_component(self, compspec,frac):
         '''Adds an element to the core, setting its fraction to frac. Adding multiple elements\n
         will change this fraction.'''
         for i in range(self.base.size):
            for j in range(self.base.size):
                if self.signal:
                     self.core.data[i][j]=self.core.data[i][j]*(1-frac)+compspec.data*frac*self.base.core[i,j];
                else:
                     self.core[i][j]=self.core[i][j]*(1-frac)+compspec.data*frac*self.base.core[i,j];
                     
     def add_shell_component(self,compspec,frac):
         '''Adds an element to the shell, setting its fraction to frac. Adding multiple elements\n
         will change this fraction.'''
         for i in range(self.base.size):
             for j in range(self.base.size):
                 if self.signal:
                     self.shell.data[i][j]=self.shell.data[i][j]*(1-frac)+compspec.data*frac*self.base.shell[i,j];
                 else:
                     self.shell[i][j]=self.shell[i][j]*(1-frac)+compspec.data*frac*self.base.shell[i,j];

    
class CoreShellBack:
    def __init__(self,a,spec,dens,signal=True):
        self.signal=signal;
        spec=spec.dens
        l=a.size()
        L=len(spec)
        mat=np.zeros((l,l))
        if signal:
            for i in range(l):
                for j in range(l):
                    mat[i][j][0:L]=hs.signals.Signal1D(spec);
        else:
            for i in range(l):
                for j in range(l):
                    mat[i][j][0:L]=spec;            
        self.size=a.size;
        self.base=a.base;
        self.back=mat
        
    def get_backg(self):
        return self.back;
    def is_signal(self):
        return self.signal

#%% sim/coreshellp_3D.py
"""
Created on Tue Feb  9 14:20:11 2021

@author: Jonas
"""

def cut_spectrum_bottom(hz,at):
    at=float(at);
    hz.isig[0:at].data-=hz.isig[0:at].data;
    return None
def cut_spectrum_top(hz,at):
    at=float(at);
    hz.isig[at:-1].data-=hz.isig[at:-1].data;
    return None
def cut_spectrum_range(hz,start,stop):
    start,stop=float(start),float(stop);
    hz.isig[start:stop].data-=hz.isig[start:stop];
    return None

class EdxMat_3D:
    def __init__(self,size,r: float,dens:float,l: float):
        '''size: defines the 'picture' size as size x size pixels\n
        r: defines the radius of the spherical particle\n
        dens: density of the material the particle is made of\n
        l: pixel size in the same unit as r (i.e. the area is l^2)'''
        self.size=size;
        self.r=r;
        self.dens=dens;
        self.l=l;
        rred=r/l;
        self.mid= float((size-1)/2);#-1 because indices start from zero
        #The constructor now constructs a matrix using the input arguments
        mat=np.zeros((size,size,size));
        self.mat=mat;
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    if rred**2>=((i-self.mid)**2+(j-self.mid)**2+(k-self.mid)**2):
                        mat[i][j][k]=1
        
    def __thick(self,n,m,p):
        nr=n-self.mid;
        mr=m-self.mid;
        pr=p-self.mid;
        l=self.l;
        r=self.r;
        return math.sqrt(r**2-(nr)**2-(mr)**2)*l**2; 
    

class CoreShellP_3D:
    def __init__(self,size,r1:float,r2:float,dens1:float,dens2:float,l:float):
        '''Constructs an object modeling a core-shell spherical particle.\n
        The object properties core and shell  are matrices \n
        modeling the core and shell respectively.\n
        All arguments are defined as for EdxMat.'''
        
        x=EdxMat_3D(size,r1,dens1,l);
        y=EdxMat_3D(size,r2,dens2,l);
        self.size=size;
        self.l=l;
        if r1>=r2:
            z=EdxMat_3D(size,r2,dens1,l)
            self.shell=x.mat-z.mat;
            self.core=y.mat;

        else:
            z=EdxMat_3D(size,r1,dens2,l)
            self.shell=y.mat-z.mat;
            self.core=x.mat;


class CoreShellSpec_3D:
     def __init__(self,a,spec1,spec2,signal=True):
         '''Takes a CoreShellP object a and inserts the core spectrum \n
         spec1 and the shell spectrum spec2 at the appropriate places.\n
         Both spectra must be HyperSpy 1D signals.'''
         self.signal=signal
         L=len(spec1.data)
         arr=np.empty((a.size,a.size,a.size,L))
         core=arr.copy()
         shell = arr.copy()
         for i in range(0,a.size):
             for j in range(0,a.size):
                 for k in range(0,a.size):
                     core[i,j,k,0:L]=a.core[i,j,k]*spec1.data
                     shell[i,j,k,0:L]=a.shell[i,j,k]*spec2.data
         if signal:
             self.core = hs.signals.Signal1D(core)
             self.shell = hs.signals.Signal1D(shell)
         else:
             self.core=core;
             self.shell=shell;
         self.full=self.core+self.shell   #Ensures backwards compatibility
         self.base=a;
     def getmatr(self):
         '''Returns the matrix representing the entire core-shell particle'''
         return self.core+self.shell
     def is_signal(self):
         return self.signal;
     def add_background(self,backspec, thickness):
         '''Adds a background to the core-shell particle. The thickness is a multiplier\n
         on the original intensity.'''
         for i in range(self.base.size):
             for j in range(self.base.size):
                 for k in range(self.base.size):
                     if self.signal:
                         self.core.data[i][j][k]+=backspec.data*thickness;
                     else:
                         self.core[i][j][k]+=backspec.data*thickness;
                     
     def add_core_component(self, compspec,frac):
         '''Adds an element to the core, setting its fraction to frac. Adding multiple elements\n
         will change this fraction.'''
         for i in range(self.base.size):
            for j in range(self.base.size):
                if self.signal:
                     self.core.data[i][j]=self.core.data[i][j]*(1-frac)+compspec.data*frac*self.base.core[i,j];
                else:
                     self.core[i][j]=self.core[i][j]*(1-frac)+compspec.data*frac*self.base.core[i,j];
                     
     def add_shell_component(self,compspec,frac):
         '''Adds an element to the shell, setting its fraction to frac. Adding multiple elements\n
         will change this fraction.'''
         for i in range(self.base.size):
             for j in range(self.base.size):
                 if self.signal:
                     self.shell.data[i][j]=self.shell.data[i][j]*(1-frac)+compspec.data*frac*self.base.shell[i,j];
                 else:
                     self.shell[i][j]=self.shell[i][j]*(1-frac)+compspec.data*frac*self.base.shell[i,j];

    
