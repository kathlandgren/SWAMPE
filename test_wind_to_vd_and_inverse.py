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
a1=np.pi/4#p.a1
test=1
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
# K4=K4*10**10
dt=100 #dt/10

# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

SU0, sina, cosa, etaamp,Phiamp=ic.test1_init(a, omega, a1)
etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,test,etaamp,a,sina,cosa,Phibar,Phiamp)
Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)


fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fmn[0,1]=omega/np.sqrt(0.375)

flatlon=np.zeros((J,I))
for i in range(I):
    flatlon[:,i]=2*omega*mus

tstepcoeffmn=tstep.tstepcoeffmn(M,N,a)
tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,a)
tstepcoeff2=tstep.tstepcoeff2(J,M,dt,a)
mJarray=tstep.mJarray(J,M)
marray=tstep.marray(M, N)
narray=tstep.narray(M,N)


def wind_test(U,V,I,J,M,N,Pmn,Hmn,w,tstepcoeff,tstepcoeffmn,mJarray,marray,dt,fmn):
   ## This test takes a wind field, converts it to vorticity and divergence and recreates the wind field 
    Um=rfl.fwd_fft_trunc(U,I,M)
    Vm=rfl.fwd_fft_trunc(V,I, M)

    eta,delta,etamn,deltamn=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)    #why does it need a dt?
        
    Unew,Vnew=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
        #plotting
    plt.contourf(lambdas, mus, eta)
    plt.colorbar()
    plt.title('eta')
    plt.show()    
        
        
    plt.contourf(lambdas, mus, U-Unew)
    plt.colorbar()
    plt.title('U error')
    plt.show()
    
    #plotting
    plt.contourf(lambdas, mus, U)
    plt.colorbar()
    plt.title('U IC')
    plt.show()
    
    plt.contourf(lambdas, mus, Unew)
    plt.colorbar()
    plt.title('U Transform')
    plt.show()
    
    plt.contourf(lambdas, mus, V-Vnew)
    plt.colorbar()
    plt.title('V error')
    plt.show()
    
    #plotting
    plt.contourf(lambdas, mus, V)
    plt.colorbar()
    plt.title('V IC')
    plt.show()
    
    plt.contourf(lambdas, mus, Vnew)
    plt.colorbar()
    plt.title('V Transform')
    plt.show()
    
    return Unew, Vnew


## This test takes vorticity and divergence, converts them to a wind field and recreates vorticity and divergence


def vor_div_test(eta,delta,I,J,M,N,Pmn,Hmn,w,tstepcoeff,tstepcoeffmn,mJarray,marray,dt,fmn):

    # Um=rfl.fwd_fft_trunc(U,I,M)
    # Vm=rfl.fwd_fft_trunc(V,I, M)

    # eta,delta,etamn,deltamn=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)    #why does it need a dt?
    
    deltam=rfl.fwd_fft_trunc(deltaic0,I,M)
    deltamn=rfl.fwd_leg(deltam,J,M,N,Pmn,w)

    etam=rfl.fwd_fft_trunc(etaic0,I,M)
    etamn=rfl.fwd_leg(etam,J,M,N,Pmn,w)

    U,V=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
    Um=rfl.fwd_fft_trunc(U,I,M)
    Vm=rfl.fwd_fft_trunc(V,I, M)
    
    etanew,deltanew,etamnnew,deltamnnew=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)    #why does it need a dt?
    
        #plotting
    plt.contourf(lambdas, mus, eta)
    plt.colorbar()
    plt.title('eta IC')
    plt.show()    
        
        
    plt.contourf(lambdas, mus,etanew)
    plt.colorbar()
    plt.title('eta new')
    plt.show()
    
    #plotting
    plt.contourf(lambdas, mus, eta-etanew)
    plt.colorbar()
    plt.title('eta error')
    plt.show()
    
    return etanew, deltanew, etamnnew, deltamnnew,U,V

## This is a test for the Laplacian

## This is a test from Liu and Showman 
def inverse_wind_test(U,V,etaic0,deltaic0,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray):
    #Tests eq. 5.24-5.25 in Hack and Jakob, trying to get the initial winds from the initial vorticity and divergence.
    deltam=rfl.fwd_fft_trunc(deltaic0,I,M)
    deltamn=rfl.fwd_leg(deltam,J,M,N,Pmn,w)
    # deltamntemp=deltamn
    # deltamn[0:M,:]=deltamntemp[1:M+1,:]
    deltamn[:,0]=0

    etam=rfl.fwd_fft_trunc(etaic0,I,M)
    etamn=rfl.fwd_leg(etam,J,M,N,Pmn,w)
    

    # etamntemp=etamn
    # etamn[0:M,:]=etamntemp[1:M+1,:]
    
    etamn[:,0]=0
    #etamn[0,1]=fmn[0,1]
    # etamnBig=np.zeros((M+2,M+2))
    # etamnBig[0:M+1,1:M+2]=etamn
    
    Um=rfl.fwd_fft_trunc(U,I,M)
    Vm=rfl.fwd_fft_trunc(V,I, M)

    etaD,deltaD,etamnD,deltamnD=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt) 
    
    Unew,Vnew=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
        
    #plotting
    plt.contourf(lambdas, mus, U-Unew)
    plt.colorbar()
    plt.title('U error')
    plt.show()
    
    #plotting
    plt.contourf(lambdas, mus, U)
    plt.colorbar()
    plt.title('U IC')
    plt.show()
    
    plt.contourf(lambdas, mus, Unew)
    plt.colorbar()
    plt.title('U Transform')
    plt.show()
    
        
    plt.contourf(np.arange(M+1), np.arange(M+1), deltamnD)
    plt.colorbar()
    plt.title('deltamn diagnostic')
    plt.show()
    
    plt.contourf(np.arange(M+1), np.arange(M+1), deltamn)
    plt.colorbar()
    plt.title('deltamn original')
    plt.show()
    
    plt.contourf(np.arange(M+1), np.arange(M+1), deltamn-deltamnD)
    plt.colorbar()
    plt.title('deltamn difference')
    plt.show()
    
    
    # plt.contourf(lambdas, mus, V-Vnew)
    # plt.colorbar()
    # plt.title('V error')
    # plt.show()
    
    # #plotting
    # plt.contourf(lambdas, mus, V)
    # plt.colorbar()
    # plt.title('V IC')
    # plt.show()
    
    # plt.contourf(lambdas, mus, Vnew)
    # plt.colorbar()
    # plt.title('V Transform')
    # plt.show()
    
    return Unew, Vnew, etamnD, etamn

def diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt):
    coeff=tstepcoeff/(2*dt)
    etacomp1prep=np.multiply(np.multiply(coeff,(1j)*mJarray),Vm)
    etacomp2prep=np.multiply(coeff,Um)
    
    zetamn=rfl.fwd_leg(etacomp1prep, J, M, N, Pmn, w)+rfl.fwd_leg(etacomp2prep, J, M, N, Hmn, w)
    
    etamn=zetamn+fmn
    
    deltacomp1prep=np.multiply(np.multiply(coeff,(1j)*mJarray),Um)
    deltacomp2prep=np.multiply(coeff,Vm)
    
    deltacomp1=rfl.fwd_leg(deltacomp1prep, J, M, N, Pmn, w)
    deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Hmn, w)
    
    deltamn=deltacomp1-deltacomp2
    
    test,newdeltam=rfl.invrs_leg(deltamn, I,J, M, N, Pmn)
    newdelta=rfl.invrs_fft(newdeltam, I)
    
    test,newetam=rfl.invrs_leg(etamn, I,J, M, N, Pmn)
    neweta=rfl.invrs_fft(newetam, I)
   
    return neweta,newdelta,etamn,deltamn

#Unew,Vnew=wind_test(Uic,Vic,I,J,M,N,Pmn,Hmn,w,tstepcoeff,tstepcoeffmn,mJarray,marray,dt,fmn)


#Unew,Vnew,etamnD, etamn=inverse_wind_test(Uic,Vic,etaic0,deltaic0,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)

etanew, deltanew, etamnnew, deltamnnew,Unew,Vnew=vor_div_test(etaic0,deltaic0,I,J,M,N,Pmn,Hmn,w,tstepcoeff,tstepcoeffmn,mJarray,marray,dt,fmn)
           
V=Vic
U=Uic

#plotting
plt.contourf(lambdas, mus, U-Unew)
plt.colorbar()
plt.title('U error')
plt.show()

#plotting
plt.contourf(lambdas, mus, U)
plt.colorbar()
plt.title('U IC')
plt.show()

plt.contourf(lambdas, mus, Unew)
plt.colorbar()
plt.title('U Transform')
plt.show()

plt.contourf(lambdas, mus, V-Vnew)
plt.colorbar()
plt.title('V error')
plt.show()

#plotting
plt.contourf(lambdas, mus, V)
plt.colorbar()
plt.title('V IC')
plt.show()

plt.contourf(lambdas, mus, Vnew)
plt.colorbar()
plt.title('V Transform')
plt.show()
#winds to vorticity and div and compare