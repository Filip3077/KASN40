# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:10:55 2021

@author: Jonas
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec
from specerr import *
from specMapDiff import *
import numpy as np
from coreshellFunctions import checkLoadFit
#%%Simuleringsresulat från 02032021
carray=np.array([[1.38072285, 1.39493912, 1.39093662, 0.63703039, 0.5811263 ,
        0.53084397, 0.520143  , 0.5129397 , 0.49666035, 0.49626743,
        0.46368683],
       [1.3998878 , 1.36791138, 1.37277823, 1.38947435, 1.42885308,
        0.57801761, 0.50321906, 0.51735107, 0.48368095, 0.45198326,
        1.95034447],
       [1.44614216, 1.36250459, 1.33988858, 1.37268687, 1.37771284,
        1.42567427, 0.51493073, 0.49342208, 0.46254094, 1.68717152,
        1.80603422],
       [1.38029739, 1.36920759, 1.35074102, 1.35307981, 1.37849661,
        1.38427533, 1.42451034, 0.49385249, 1.50615046, 1.58154548,
        1.70231599],
       [1.43453271, 1.39383058, 1.34590875, 1.33418208, 1.35772428,
        1.35526881, 1.39550051, 1.43492165, 1.51582858, 1.49670907,
        1.61637768],
       [1.48130086, 1.42886732, 1.3582567 , 1.32520328, 1.33343124,
        1.35966157, 1.37466702, 1.41106065, 1.45651039, 1.56250259,
        1.53159672],
       [1.52538137, 1.41892408, 0.49856826, 1.32939995, 1.32911108,
        1.33517297, 1.35229493, 1.37970209, 1.42209117, 1.49621963,
        1.45301763],
       [1.58899532, 0.42112679, 0.49058557, 0.53792686, 1.31627984,
        1.35042036, 1.34921399, 1.34773324, 1.39410311, 1.4402318 ,
        1.39076342],
       [0.37061971, 0.45050445, 0.4431634 , 0.51852939, 1.32551206,
        1.31543128, 1.32523656, 1.36324893, 1.35576121, 1.39604359,
        1.48697075],
       [0.39279421, 0.43398649, 0.45154156, 0.50726328, 0.56077537,
        1.30544328, 1.32101084, 1.32636987, 1.35158558, 1.3649271 ,
        1.41545551],
       [0.39652467, 0.43328009, 0.46230057, 0.49022921, 0.52993694,
        0.58791246, 0.6174209 , 1.32243349, 1.33597893, 1.36744366,
        1.36823181]]);
sarray=np.array([[0.81146643, 0.83252329, 0.84869895, 1.03069733, 1.08902899,
        1.14528856, 1.19151071, 1.2206175 , 1.28376623, 1.32410647,
        1.39191379],
       [0.83547015, 0.81499705, 0.82172269, 0.83519659, 0.86901605,
        1.02517712, 1.06950218, 1.12590835, 1.17050094, 1.22772125,
        0.28800907],
       [0.87073379, 0.82966376, 0.80175803, 0.81745346, 0.82703807,
        0.85965439, 0.99240055, 1.04190884, 1.09237464, 0.32042993,
        0.27387512],
       [0.28815966, 0.88346981, 0.83839841, 0.80712679, 0.82017098,
        0.82366812, 0.84311928, 0.98400025, 0.31430402, 0.31789925,
        0.2747053 ],
       [0.29329906, 0.90868545, 0.87078949, 0.83058726, 0.82699835,
        0.80774756, 0.82204228, 0.83909887, 0.86894839, 0.32082256,
        0.28136077],
       [0.29795106, 0.9643918 , 0.90777142, 0.86345538, 0.82481916,
        0.81445145, 0.81500241, 0.8156697 , 0.83533627, 0.85525064,
        0.27523975],
       [0.28835045, 0.36211752, 1.04920298, 0.91086494, 0.86382516,
        0.83191651, 0.80172107, 0.80613661, 0.8304403 , 0.83389969,
        0.27011233],
       [0.31603943, 1.18510679, 1.10050033, 1.04959729, 0.91171798,
        0.87298562, 0.8244691 , 0.79826682, 0.80335694, 0.81815146,
        0.28093621],
       [1.3435108 , 1.25381346, 1.18088545, 1.10756076, 0.95967065,
        0.904636  , 0.8592354 , 0.82824283, 0.80642552, 0.80679083,
        0.82379302],
       [1.43185695, 1.33107357, 1.26368224, 1.18723029, 1.12568372,
        0.9399597 , 0.89924056, 0.86088029, 0.82269667, 0.80258895,
        0.80204133],
       [1.53153444, 1.44968721, 1.36860207, 1.30020007, 1.2430561 ,
        1.1701622 , 1.10517983, 0.90628434, 0.86093023, 0.83076394,
        0.79918524]]);
earray=np.array([[0.81146643, 0.83252329, 0.84869895, 1.03069733, 1.08902899,
        1.14528856, 1.19151071, 1.2206175 , 1.28376623, 1.32410647,
        1.39191379],
       [0.83547015, 0.81499705, 0.82172269, 0.83519659, 0.86901605,
        1.02517712, 1.06950218, 1.12590835, 1.17050094, 1.22772125,
        0.28800907],
       [0.87073379, 0.82966376, 0.80175803, 0.81745346, 0.82703807,
        0.85965439, 0.99240055, 1.04190884, 1.09237464, 0.32042993,
        0.27387512],
       [0.28815966, 0.88346981, 0.83839841, 0.80712679, 0.82017098,
        0.82366812, 0.84311928, 0.98400025, 0.31430402, 0.31789925,
        0.2747053 ],
       [0.29329906, 0.90868545, 0.87078949, 0.83058726, 0.82699835,
        0.80774756, 0.82204228, 0.83909887, 0.86894839, 0.32082256,
        0.28136077],
       [0.29795106, 0.9643918 , 0.90777142, 0.86345538, 0.82481916,
        0.81445145, 0.81500241, 0.8156697 , 0.83533627, 0.85525064,
        0.27523975],
       [0.28835045, 0.36211752, 1.04920298, 0.91086494, 0.86382516,
        0.83191651, 0.80172107, 0.80613661, 0.8304403 , 0.83389969,
        0.27011233],
       [0.31603943, 1.18510679, 1.10050033, 1.04959729, 0.91171798,
        0.87298562, 0.8244691 , 0.79826682, 0.80335694, 0.81815146,
        0.28093621],
       [1.3435108 , 1.25381346, 1.18088545, 1.10756076, 0.95967065,
        0.904636  , 0.8592354 , 0.82824283, 0.80642552, 0.80679083,
        0.82379302],
       [1.43185695, 1.33107357, 1.26368224, 1.18723029, 1.12568372,
        0.9399597 , 0.89924056, 0.86088029, 0.82269667, 0.80258895,
        0.80204133],
       [1.53153444, 1.44968721, 1.36860207, 1.30020007, 1.2430561 ,
        1.1701622 , 1.10517983, 0.90628434, 0.86093023, 0.83076394,
        0.79918524]]);
corrat=np.linspace(0,1,11);
ratios=np.linspace(0,1,11);
#%%Test av 3D-plot
K,I=np.meshgrid(corrat,ratios);
fig1 = plt.figure(100)
ax = fig1.gca(projection='3d')
surf = ax.plot_surface(K, I, carray, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)
ax.set_zlim(0, None)
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
#%%Successiva sammansättningsplottar
for i in range(len(ratios)):
    plt.figure(i)
    plt.plot(ratios, carray[i])
    nr=str(10*i)
    plt.xlabel('Fraction Cu in shell')
    plt.ylabel('Relative error in core')
    plt.title(nr+'% Cu in core')
    plt.ylim(0,2)