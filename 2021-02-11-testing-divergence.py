# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:46:42 2021

@author: ek672
"""

## Import statements
# Import python packages
import numpy as np
import matplotlib.pyplot as plt
#import math

# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
#import forcing
import filters
import schwartztrauber as S


#import parameters

M = p.M
# N = p.N
# Sine of the latitude
# mus = p.mus
# Wieghts for integrating
# w = p.w
# lambdas=p.lambdas
g=p.g
taurad=p.taurad
taudrag=p.taudrag
Phibar=p.Phibar
Dheq=p.Dheq
omega=p.omega
a=p.a
a1=p.a1
test=p.test

#normalization for the spherical harmonics
normnum = 1

#colorbar settings for plotting
minlevel=p.minlevel
maxlevel=p.maxlevel



N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT=S.mnMATgen(I,J,M,N,mus)


#create a strictly zonal flow in lat-lon



 #initialize the velocity array
U=np.full((J,I),10.0)
V=np.full((J,I),0.0)

truediv=np.full((J,I),0.0)

for i in range(I):
    for j in range(J):
        U[j,i]=np.sin(lambdas[i])
        V[j,i]=0#2*lambdas[i]*mus[j]
        
        truediv[j,i]=2*np.cos(lambdas[i])/a

#generate the divergence coefficients
brmn,bimn=S.A22_A23(U,V,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)

divergence= S.A15(brmn,bimn,nMAT1,mus,M,J,normnum)/a
#divergence=np.multiply(div1,musMAT)
print(divergence-truediv)
error=divergence-truediv

plt.contourf(lambdas, mus, error)
plt.colorbar()
plt.title('error')
plt.show()

plt.contourf(lambdas, mus, truediv)
plt.colorbar()
plt.title('true div')
plt.show()

plt.contourf(lambdas, mus, divergence)
plt.colorbar()
plt.title('our divergence')
plt.show()


