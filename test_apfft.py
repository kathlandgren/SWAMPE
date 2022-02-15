# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 09:59:01 2021

@author: ek672
"""

import numpy as np
import mpmath as mp
import apfft.fft

import fft_legendre_trans as rfl


mp.dps=50

M=42
J=64
I=128

z=mp.zeros(J,I)

z[0,0]=mp.pi

np.fft.fft(z)

# error: cannot copy sequence with size 64 to array axis with dimension 8192



