# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:15:02 2021

@author: ek672
"""
import numpy as np
import matplotlib.pyplot as plt

import spectral_transform as rfl
import params as p
import initial_conditions as ic
import time_stepping as tstep


def init_test():
    #tests N, I, J initialization from M
    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(42)
    
    assert N==42, "N should be 42"
    assert I==128, "I should be 128"
    assert J==64, "J should be 64"   
    
def Pmn_Hmn_test():
    #tests Legendre polynomials against analytical definition
    M=42
    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
    Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
    
    Pmncheck=0.25*np.sqrt(15)*(1-mus**2)
    Hmncheck=0.5*np.sqrt(6)*(1-mus**2)
    
    assert np.allclose(Pmn[:,2,2],Pmncheck,atol=1e-12),"Associated Legendre polynomial for m=2, n=2 not correct"
    assert np.allclose(Hmn[:,0,1],Hmncheck,atol=1e-12),"Associated Legendre polynomial derivative for m=0, n=1 not correct"

def spectral_transform_test():
    #tests the forward spectral transform on Coriolis force
    M=106
    omega=3.2*10**(-5)
    N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
    Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
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
    a1=np.pi/2
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
    M = 106

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
    Pmn_Hmn_test()
    spectral_transform_test()
    spectral_transform_forward_inverse_test()
    wind_transform_test()
    vorticity_divergence_transform_test()
    print('All tests passed!')
    