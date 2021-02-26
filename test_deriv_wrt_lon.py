# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:17:48 2021

@author: ek672

This file tests the derivative with respect to longitude
"""

#import and initialize 

import numpy as np
import matplotlib.pyplot as plt

import fft_legendre_trans as rfl
import params as p
import initial_conditions as ic
import tstepping_new as tstep

M = 63
# N = p.N
# Length of the run in time steps
tmax = p.tmax
# Sine of the latitude
# mus = p.mus
# Wieghts for integrating
# w = p.w
# lambdas=p.lambdas
g=p.g
taurad=p.taurad
taudrag=p.taudrag
Phibar=p.Phibar
omega=p.omega
a=p.a
a1=p.a1
test=p.test


N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
# K4=K4*10**10
dt=100 #dt/10

# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

SU0, sina, cosa, etaamp,Phiamp=ic.test1_init(a, omega, a1)
etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,a,sina,cosa,etaamp,test)
Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)


fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fmn[0,1]=omega/np.sqrt(0.375)

tstepcoeffmn=tstep.tstepcoeffmn(M,N,a)
tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,a)
tstepcoeff2=tstep.tstepcoeff2(J,M,dt,a)
mJarray=tstep.mJarray(J,M)
marray=tstep.marray(M, N)
narray=tstep.narray(M,N)


#Make a J by I vector field

testfield=np.zeros([J,I])
for i in range(I):
    testfield[:,i]=np.sin(lambdas[i])
    
plt.contourf(lambdas, mus, testfield)
plt.colorbar()
plt.title('Before derivation')
plt.show()
    

#Make its derivative wrt I
#divide by a(1-mu^2)
testfieldderiv=np.zeros([J,I])
for i in range(I):
    for j in range(J):
        testfieldderiv[j,i]=np.cos(lambdas[i])*(1/(a*(1-mus[j]**2)))    
plt.contourf(lambdas, mus, testfieldderiv)
plt.colorbar()
plt.title('True derivative')
plt.show()


#take the derivative in spectral space
testm=rfl.fwd_fft_trunc(testfield,I,M)
coeff=tstepcoeff/(2*dt) #how it's done in the code -- probably error is introduced
derivprep=np.multiply(np.multiply(coeff,(1j)*mJarray),testm)

derivmn=rfl.fwd_leg(derivprep, J, M, N, Pmn, w)

#invert the derivative

test,derivm=rfl.invrs_leg(derivmn, I,J, M, N, Pmn)
derivnew=rfl.invrs_fft(derivm, I)

#plot the difference between derivatives
plt.contourf(lambdas, mus, derivnew)
plt.colorbar()
plt.title('New derivative')
plt.show()

plt.contourf(lambdas, mus, testfieldderiv-derivnew)
plt.colorbar()
plt.title('True-new error')
plt.show()
