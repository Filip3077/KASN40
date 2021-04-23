# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 08:19:53 2021

@author: Jonas
"""

import matplotlib.pyplot as plt
import hyperspy.api as hs
from coreshellp import *
from specerr import *
from specMapDiff import *
import numpy as np
from loadassign import checkLoadFit
from scipy.ndimage import gaussian_filter
from k_factors import silver_k_factor
from radial_profile import transfer_elements

p22_1=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 22 nm\\EDS Data_particle 1.rpl',signal_type='EDS_TEM').T
p22_2=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 22 nm\\EDS Data_particle 2.rpl',signal_type='EDS_TEM').T
p22_3=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 22 nm\\EDS Data_particle 3.rpl',signal_type='EDS_TEM').T
p22=[p22_1,p22_2,p22_3];
p30_1=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 1.rpl',signal_type='EDS_TEM').T
p30_2=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 2.rpl',signal_type='EDS_TEM').T
p30_3=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 3.rpl',signal_type='EDS_TEM').T
p30_4=hs.load('C:\\Users\\Jonas\\Documents\\KASN40 Project course\\Experimentell data\\AuZn 30 nm\\EDS Data_particle 4.rpl',signal_type='EDS_TEM').T
p30=[p30_1,p30_2,p30_3,p30_4];
#%%Plotta
p22[1].decomposition('sklearn_pca')
p22[1].plot_explained_variance_ratio
    



