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

M = 42
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
    #tests N, I, J initialization from M
    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(42)
    
    assert N==42, "N should be 42"
    assert I==128, "I should be 128"
    assert J==64, "J should be 64"   
    
#def Pmn_test():

def spectral_transform_test():
    #tests the forward spectral transform on Coriolis force
    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
    f=np.zeros((J,I))
    for i in range(I):
             for j in range(J):
                 f[j,i]=2*omega*mus[j]

    #fmn=np.zeros([M+1,N+1]) 
    fm=np.zeros((J,M+1))
    fm=rfl.fwd_fft_trunc(f, I, M)

    fmn=rfl.fwd_leg(fm, J, M, N, Pmn, w)
    
    fmncheck=np.zeros([M+1,N+1]) 
    fmncheck[0,1]=omega/np.sqrt(0.375)
    
    assert np.allclose(fmn, fmncheck, atol=1e-12), "Spectral coefficients should be close"
   



def spectral_transform_forward_inverse_test():
    ## This test takes a scalar field, transforms is to spectral space and back
    
    M = 63

    
    omega=p.omega
    a=p.a
    a1=np.pi/4#p.a1
    test=1


    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

    # Associated Legendre Polynomials and their derivatives
    Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

    SU0, sina, cosa, etaamp,Phiamp=ic.test1_init(a, omega, a1)
    etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,a,sina,cosa,etaamp,test)
    Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
    #Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)
    
    Uicm=rfl.fwd_fft_trunc(Uic,I,M)

     
    Uicmn=rfl.fwd_leg(Uicm,J,M,N,Pmn,w)
        
    Uicmnew=rfl.invrs_leg(Uicmn,I,J,M,N,Pmn)
    
    Uicnew=rfl.invrs_fft(Uicmnew,I)
    
    assert np.allclose(Uic, Uicnew, atol=10**(-11)), "Error too high"



def wind_transform_test():
   ## This test takes a wind field, converts it to vorticity and divergence and recreates the wind field 
    M = 63

    omega=p.omega
    a=p.a
    a1=np.pi/4#p.a1
    test=1
    dt=30

    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

    # Associated Legendre Polynomials and their derivatives
    Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
    
    fmn=np.zeros([M+1,N+1]) 
    fmn[0,1]=omega/np.sqrt(0.375)
    
    tstepcoeffmn=tstep.tstepcoeffmn(M,N,a)
    tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,a)
    mJarray=tstep.mJarray(J,M)
    marray=tstep.marray(M, N)


    SU0, sina, cosa, etaamp,Phiamp=ic.test1_init(a, omega, a1)
    etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,a,sina,cosa,etaamp,test)
    U,V=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
   
    Um=rfl.fwd_fft_trunc(U,I,M)
    Vm=rfl.fwd_fft_trunc(V,I, M)
 
    eta,delta,etamn,deltamn=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)
        
    Unew,Vnew=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
    assert np.allclose(U, Unew, atol=10**(-11)), "U error too high"
    assert np.allclose(V, Vnew, atol=10**(-11)), "V error too high"

def vorticity_divergence_transform_test():

## This test takes vorticity and divergence, converts them to a wind field and recreates vorticity and divergence
    
    M = 63

    omega=p.omega
    a=p.a
    a1=np.pi/4#p.a1
    test=1
    dt=30

    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

    # Associated Legendre Polynomials and their derivatives
    Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
    
    fmn=np.zeros([M+1,N+1]) 
    fmn[0,1]=omega/np.sqrt(0.375)
    
    tstepcoeffmn=tstep.tstepcoeffmn(M,N,a)
    tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,a)
    mJarray=tstep.mJarray(J,M)
    marray=tstep.marray(M, N)


    SU0, sina, cosa, etaamp,Phiamp=ic.test1_init(a, omega, a1)
    etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,a,sina,cosa,etaamp,test)
    U,V=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
    deltam=rfl.fwd_fft_trunc(deltaic0,I,M)
    deltamn=rfl.fwd_leg(deltam,J,M,N,Pmn,w)

    etam=rfl.fwd_fft_trunc(etaic0,I,M)
    etamn=rfl.fwd_leg(etam,J,M,N,Pmn,w)

    U,V=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
    Um=rfl.fwd_fft_trunc(U,I,M)
    Vm=rfl.fwd_fft_trunc(V,I, M)
    
    etanew,deltanew,etamnnew,deltamnnew=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)
    
    assert np.allclose(etaic0, etanew, atol=10**(-11)), "Vorticity eta error too high"
    assert np.allclose(deltaic0, deltanew, atol=10**(-11)), "Divergence delta error too high"
    

if __name__=='__main__':
    init_test()
    spectral_transform_test()
    spectral_transform_forward_inverse_test()
    wind_transform_test()
    vorticity_divergence_transform_test()
    print('All tests passed')
    


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