# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 21:39:32 2020

@author: Martin
"""

import hyperspy.api as hs
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

s = hs.load('ROI4_2byte.rpl').as_signal1D(0)

s.set_signal_type('EDS_TEM')
s.axes_manager[-1].name = 'E'
s.axes_manager['E'].scale = 0.01
s.axes_manager['E'].offset = -0.2
s.axes_manager.signal_axes[0].units = 'keV' 

keV = np.arange(-0.2, 13.81, 0.01)

def varimax(Phi , gamma = 1.0 , q = 100  , tol = 1e-6):
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
        
        if d/d_old < (1+tol) :
            break
    return R


#%% Cut and rebin
s = s.inav[44:512,178:230] #43:512 and 177:230 original

bin_scale = [2, 52, 1]

s1 = s.rebin(scale=bin_scale)
s1.change_dtype('float64')

#%% PCA part

s1.decomposition(normalize_poissonian_noise=True)
s1.plot_explained_variance_ratio(n=20, vline=True)

sc1 = s1.get_decomposition_model(4)
res1 = (s1 - sc1)

res_map=res1.sum(-1).data[0,:]
res_spec=res1.sum([0,1]).data

s1_factors = s1.get_decomposition_factors()
s1_loadings = s1.get_decomposition_loadings().sum(-1)

# s1_factors.inav[0].isig[0:180].plot() plot low energy region of factor

#%% Varimax part

nfac = 4

s1_factors_selected = s1_factors.data[0:nfac,:]
s1_loadings_selected = s1_loadings.data[0:nfac,:]

R =varimax(s1_loadings_selected.T)

s1_loadings_selected_rot = np.matmul(s1_loadings_selected.T, R).T
s1_factors_selected_rot = np.matmul(linalg.inv(R), s1_factors_selected)

s1_factors_selected_rot[2,:] = -s1_factors_selected_rot[2,:]
s1_loadings_selected_rot[2,:] = -s1_loadings_selected_rot[2,:]

#s1_factors_selected_rot[4,:] = -s1_factors_selected_rot[4,:]
#s1_loadings_selected_rot[4,:] = -s1_loadings_selected_rot[4,:]

#%% Plotting

fig1, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))

ax1.plot(keV, s1_factors_selected_rot[0,:], label='C1')
ax1.plot(keV, s1_factors_selected_rot[1,:], label='C2')
ax1.plot(keV, s1_factors_selected_rot[2,:], label='C3')
ax1.plot(keV, s1_factors_selected_rot[3,:], label='C4')
#ax1.plot(keV, s1_factors_selected_rot[4,:], label='C5')
ax1.axis([-0.2,1.6,0, 10])
ax1.legend()

ax2.plot(keV, s1_factors_selected_rot[0,:], label='C1')
ax2.plot(keV, s1_factors_selected_rot[1,:], label='C2')
ax2.plot(keV, s1_factors_selected_rot[2,:], label='C3')
ax2.plot(keV, s1_factors_selected_rot[3,:], label='C4')
#ax2.plot(keV, s1_factors_selected_rot[4,:], label='C5')
ax2.axis([1.6,12.5,0, 40])
ax2.legend()

fig2, ax3 = plt.subplots()

ax3.plot(s1_loadings_selected_rot[0,:], label='C1')
ax3.plot(s1_loadings_selected_rot[1,:], label='C2')
ax3.plot(s1_loadings_selected_rot[2,:], label='C3')
ax3.plot(s1_loadings_selected_rot[3,:], label='C4')
#ax3.plot(s1_loadings_selected_rot[4,:], label='C5')
ax3.axis([0,233,-15,55])
ax3.legend()

#%% Comparisons with WC and W2C

WC_spec = hs.load("WC.msa")
WCdata = WC_spec.isig[:1401].data
WCdata = WCdata/max(WCdata[840:880])*max(s1_factors_selected_rot[0,840:880])

W2C_spec = hs.load("W2C.msa")
W2Cdata = W2C_spec.isig[:1401].data
W2Cdata = W2Cdata/max(W2Cdata[840:880])*max(s1_factors_selected_rot[0,840:880])

fig3, ax4 = plt.subplots()
ax4.plot(keV, WCdata+0.4, label='WC sim.')
ax4.plot(keV, W2Cdata+0.4, label='W2C sim.')
ax4.plot(keV, s1_factors_selected_rot[0,:], label='C1')
ax4.plot(keV, s1_factors_selected_rot[1,:], label='C2')
ax4.axis([-0.2,1.6,0, 3 ])
ax4.legend()