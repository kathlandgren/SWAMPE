# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:15:02 2021

@author: ek672
"""
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


def transform_test(Xi,I,J,M,N,Pmn,w,mus,lambdas):
    ## This test takes a scalar field, transforms is to spectral space and back
    Xim=rfl.fwd_fft_trunc(Xi,I,M)

     
    Ximn=rfl.fwd_leg(Xim,J,M,N,Pmn,w)
        
    temp, Ximnew=rfl.invrs_leg(Ximn,I,J,M,N,Pmn)
    print(np.shape(Ximnew))
    
    Xinew=rfl.invrs_fft(Ximnew,I)
    
    #plotting
    plt.contourf(lambdas, mus, Xi-Xinew)
    plt.colorbar()
    plt.title('error')
    plt.show()
    
    #plotting
    plt.contourf(lambdas, mus, Xi)
    plt.colorbar()
    plt.title('IC')
    plt.show()
    
    plt.contourf(lambdas, mus, Xinew)
    plt.colorbar()
    plt.title('Transform')
    plt.show()
    
    return Xinew


Unew=transform_test(Uic,I,J,M,N,Pmn,w,mus,lambdas)

