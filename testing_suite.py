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
Dheq=p.DPhieq
omega=p.omega
a=p.a
a1=p.a1
test=p.test


N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
# K4=K4*10**10
dt=100 #dt/10

# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

#SU0, sina, cosa, etaamp=ic.test1_init(a, omega, a1)
#etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,a,sina,cosa,etaamp,test)
#Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
#Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)


# fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
# fmn[0,1]=omega/np.sqrt(0.375)

# tstepcoeffmn=tstep.tstepcoeffmn(M,N,a)
# tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,a)
# tstepcoeff2=tstep.tstepcoeff2(J,M,dt,a)
# mJarray=tstep.mJarray(J,M)
# marray=tstep.marray(M, N)
# narray=tstep.narray(M,N)


def init_test():
    
    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(42)
    
    assert N==42, "N should be 42"
    assert I==128, "I should be 128"
    assert J==64, "J should be 64"   


if __name__=='__main__':
    init_test()
    print('All tests passed')
    

# def transform_test(U,I,J,M,N,Pmn,w,mus,lambdas):
#     ## This test takes a scalar field, transforms is to spectral space and back
#     Um=rfl.fwd_fft_trunc(U,I,M)

     
#     Umn=rfl.fwd_leg(Um,J,M,N,Pmn,w)
        
#     temp, Umnew=rfl.invrs_leg(Umn,I,J,M,N,Pmn)
#     print(np.shape(Umnew))
    
#     Unew=rfl.invrs_fft(Umnew,I)
    
#     #plotting
#     plt.contourf(lambdas, mus, U-Unew)
#     plt.colorbar()
#     plt.title('U error')
#     plt.show()
    
#     #plotting
#     plt.contourf(lambdas, mus, U)
#     plt.colorbar()
#     plt.title('U IC')
#     plt.show()
    
#     plt.contourf(lambdas, mus, Unew)
#     plt.colorbar()
#     plt.title('Transform')
#     plt.show()
    
#     return Unew

# def wind_test(U,V,I,J,M,N,Pmn,Hmn,w,tstepcoeff,tstepcoeffmn,mJarray,marray,dt,fmn):
#    ## This test takes a wind field, converts it to vorticity and divergence and recreates the wind field 
#     Um=rfl.fwd_fft_trunc(U,I,M)
#     Vm=rfl.fwd_fft_trunc(V,I, M)

#     eta,delta,etamn,deltamn=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)    #why does it need a dt?
        
#     Unew,Vnew=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
#         #plotting
#     plt.contourf(lambdas, mus, U-Unew)
#     plt.colorbar()
#     plt.title('U error')
#     plt.show()
    
#     #plotting
#     plt.contourf(lambdas, mus, U)
#     plt.colorbar()
#     plt.title('U IC')
#     plt.show()
    
#     plt.contourf(lambdas, mus, Unew)
#     plt.colorbar()
#     plt.title('Transform')
#     plt.show()
    
#     return Unew, Vnew
## This test takes vorticity and divergence, converts them to a wind field and recreates vorticity and divergence

## This is a test for the Laplacian

## This is a test from Liu and Showman 


# Unew=transform_test(Uic,I,J,M,N,Pmn,w,mus,lambdas)

# Unew,Vnew=wind_test(Uic,Vic,I,J,M,N,Pmn,Hmn,w,tstepcoeff,tstepcoeffmn,mJarray,marray,dt,fmn)