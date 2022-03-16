# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 16:13:17 2022

@author: ek672
"""

## Import statements
# Import python packages
import numpy as np
import matplotlib.pyplot as plt

# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import forcing
import filters
import continuation as cont


import pickle

M=42
dt=180
a=p.a

N,I,J,otherdt,K4,lambdas,mus,w=ic.spectral_params(M)
K4=0.5*10**(16)
K6=1.24*10**33#1.24*10**30

#K6=10**34 #empirically the same as KEarth, I guess

Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

sigma6=filters.sigma6(M,N,K6,a, dt)
sigmaPhi6=filters.sigma6Phi(M, N, K6, a, dt)


#4th degree sigma
sigma4=filters.sigma(M,M,K4,p.a,dt)
sigma4=filters.sigmaPhi(M, M, K4, p.a, dt)
        
#Earth sigma
sigmaEarth=filters.sigma(M,M,K4,6.37122*10**(6),1200)
sigmaPhi=filters.sigmaPhi(M, M, K4, 6.37122*10**(6), 1200)


#HD4-K6 degree sigma
sigma4K6=filters.sigma(M,M,K6,p.a,dt)
sigma4K6=filters.sigmaPhi(M, M, K6, p.a, dt)

plt.plot(sigma6[0,:])
plt.plot(sigma4[0,:])
plt.plot(sigmaEarth[0,:])
# plt.plot(sigma4K6[0,:])
plt.legend(("HD6","HD4","HD4-Earth","HD4K6"))
plt.xlabel("wavenumber")
plt.ylabel("sigma factor")
#plt.ylim((0.99999,1.00001))