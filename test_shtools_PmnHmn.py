# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 12:46:59 2021

@author: ek672
"""

import numpy as np
import matplotlib.pyplot as plt

import fft_legendre_trans as rfl
import params as p
import initial_conditions as ic
import tstepping_new as tstep

M = 106
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
a1=0.05#p.a1
test=1


N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
# K4=K4*10**10
dt=100 #dt/10

# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
PmnSH, HmnSH = rfl.PmnHmnSH(J, M, N, mus)

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

m=1
n=1

plt.plot(mus,Pmn[:,m,n])
plt.title('Pmn '+str(m)+', '+str(n))
plt.show()
plt.plot(mus,Hmn[:,m,n])
plt.title('Hmn '+str(m)+', '+str(n))
plt.show()
plt.plot(mus,PmnSH[:,m,n])
plt.title('PmnSH '+str(m)+', '+str(n))
plt.show()
plt.plot(mus,HmnSH[:,m,n])
plt.title('HmnSH '+str(m)+', '+str(n))
plt.show()
plt.plot(mus,Pmn[:,m,n]-PmnSH[:,m,n])
plt.title('Pmn-PmnSH  '+str(m)+', '+str(n))
plt.show()
plt.plot(mus,Hmn[:,m,n]-HmnSH[:,m,n])
plt.title('Hmn-HmnSH  '+str(m)+', '+str(n))
plt.show()
