# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:50:02 2021

@author: Filip
"""

import math;
import numpy as np;
from edxmat import EdxMat
import hyperspy.api as hs
from coreshellp import CoreShellP, CoreShellSpec


s = hs.load("../Spectra/MC simulation of  a 0.020 µm base, 0.020 µm high block*.msa",stack=True,signal_type="EDS_TEM")

sCu = s.inav[-1]
sAg  = s.inav[0]


