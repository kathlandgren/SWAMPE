# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:06:35 2020

@author: ek672
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sp


# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import filters


import fft_legendre_trans as rfl


M=63
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

Pmn, Hmn=rfl.PmnHmnOld(J,M,N,mus)

Pmn2,Hmn2=rfl.PmnHmn(J,M,N,mus)

mtest=30
ntest=30

factor=np.max(Pmn[:,mtest,ntest])/np.max(Pmn2[:,mtest,ntest])

print(factor)

plt.plot(Hmn[:,mtest,ntest])
plt.plot(factor*Hmn2[:,mtest,ntest])
