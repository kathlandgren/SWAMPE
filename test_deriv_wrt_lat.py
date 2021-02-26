# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:17:48 2021

@author: ek672

This file tests the derivative with respect to latitude
"""

#import and initialize 

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sp

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
#Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

Pmn=np.zeros((J,M+1,N+1))
Hmn=np.zeros((J,M+1,N+1))
Hmn2=np.zeros((J,M+1,N+1))
Pmntemp=np.zeros((J,M+1,N+1))
Hmntemp=np.zeros((J,M+1,N+1))
Hmntemp2=np.zeros((J,M+1,N+1))
for j in range (0,J):
    Pmntemp[j],Hmntemp[j] = sp.lpmn(M,N,mus[j])
    Hmntemp2[j,:,:] = (1-mus[j]**2)*Hmntemp[j,:,:]
    
for m in range (0,M+1):
    for n in range (m,N+1):
        scaling_term = np.sqrt((((2*n)+1)*math.factorial(n-m))/(2*math.factorial(n+m)))
        Pmn[:,m,n] = scaling_term*Pmntemp[:,m,n]
        Hmn[:,m,n] = scaling_term*Hmntemp[:,m,n]
        Hmn2[:,m,n] = scaling_term*Hmntemp2[:,m,n]
        if m % 2 == 1:
            if n>0:
                Pmn[:,m,n]=-Pmn[:,m,n]
                Hmn[:,m,n]=-Hmn[:,m,n]
                Hmn2[:,m,n]=-Hmn2[:,m,n]

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
for j in range(J):
    testfield[j,:]=np.sin(mus[j]*np.pi/2)
    
plt.contourf(lambdas, mus, testfield)
plt.colorbar()
plt.title('Before derivation')
plt.show()
    

#Make its derivative wrt I
#divide by a(1-mu^2)
testfieldderiv=np.zeros([J,I])
for i in range(I):
    for j in range(J):
        testfieldderiv[j,i]=np.cos(mus[j])*(0.5*np.pi/a)    
plt.contourf(lambdas, mus, testfieldderiv)
plt.colorbar()
plt.title('True derivative')
plt.show()


#take the derivative in spectral space
testm=rfl.fwd_fft_trunc(testfield,I,M)
testback=rfl.invrs_fft(testm, I)
coeff=1#tstepcoeff/(2*dt) #how it's done in the code -- probably error is introduced
derivprep=np.multiply(coeff,testm)
derivprep[np.abs(derivprep)<10**(-16)]=0.0
derivmn=-rfl.fwd_leg(derivprep, J, M, N, Hmn, w)

#invert the derivative

test,derivm=rfl.invrs_leg(derivmn, I,J, M, N, Pmn)
derivnew=rfl.invrs_fft(derivm, I)
print(np.max(derivnew))
print(np.min(derivnew))

#plot the difference between derivatives
plt.contourf(lambdas, mus, derivnew)
plt.colorbar()
plt.title('New derivative')
plt.show()

plt.contourf(lambdas, mus, testfieldderiv-derivnew)
plt.colorbar()
plt.title('True-new error')
plt.show()
