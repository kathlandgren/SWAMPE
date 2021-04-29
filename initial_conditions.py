
"""
Created on Wed Jun  3 18:29:42 2020

@author: ek672

This file contains the initial conditions for the simulation
"""

import numpy as np
import scipy.special as sp

#local imports
#import params as p
#import pyshtools as pysh


def test1_init(a,omega,a1):
    #Parameters for Test 1 in Williamson et al. (1992)
    SU0=2.0*np.pi*a/(3600.0*24*12) 
    sina=np.sin(a1) #sine of the angle of advection
    cosa=np.cos(a1)  #cosine of the angle of advection
    etaamp = 2.0*((SU0/a)+omega) #relative vorticity amplitude
    Phiamp = (SU0*a*omega + 0.5*SU0**2) #geopotential height amplitude
    
    return SU0, sina, cosa, etaamp, Phiamp

def state_var_init(I,J,mus,lambdas,test,etaamp,*args):
    """Initializes the state variables
    
    :param I: number of longitudes
    :type I: int
    
    :param J: number of latitudes
    :type J: int 

    :return: J by I data arrays for eta0, eta1, delta0, delta1, and phi0, phi1
    :rtype: arrays of float64
    """
    etaic0=np.zeros((J,I))
    Phiic0=np.zeros((J,I))
    deltaic0=np.zeros((J,I))
    

    if test<=2:
        a,sina,cosa,Phibar,Phiamp=args
    if test==1:
        bumpr=a/3 #radius of the bump
        mucenter=0
        lambdacenter=3*np.pi/2
        for i in range(I):
            for j in range(J):
                etaic0[j,i]=etaamp*(-np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*sina+(mus[j])*cosa)
               
                dist=a*np.arccos(mucenter*mus[j]+np.cos(np.arcsin(mucenter))*np.cos(np.arcsin(mus[j]))*np.cos(lambdas[i]-lambdacenter))
                if dist < bumpr:
                    Phiic0[j,i]=(Phibar/2)*(1+np.cos(np.pi*dist/(bumpr)))#*p.g
    
    elif test==2: #Williamson Test 2 as documented in stswm FORTRAN implementation (see init.i)
        for i in range(I):
            for j in range(J):
                latlonarg = -np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*sina+(mus[j])*cosa
                etaic0[j,i]=etaamp*(latlonarg)

                Phiic0[j,i]=((Phibar-Phiamp)*(latlonarg)**2)#/g
    
    elif test==10 or test==11:
        for i in range(I):
            for j in range(J):
                etaic0[j,i]=etaamp*(-np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*0+(mus[j])*1)
                
                
               
    etaic1=etaic0 #need two time steps to initialize
    deltaic1=deltaic0

    Phiic1=Phiic0
    return etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1

def spectral_params(M):
    N=M
    #set dimensions according to Jakob and Hack (1993), Tables 1, 2, and 3
    if M==42:
        J=64
        I=128
        dt=1200
        K4=0.5*10**(16)
    elif M==63:
        J=96
        I=192
        dt=900
        K4=10.0**(15)
    elif M==106:
        J=160
        I=320
        dt=600
        K4=1.25*10**(14)
    elif M==170:
        J=256
        I=512
        dt=450
        K4=2.00*10*(13)
    elif M==213:
        J=320
        I=640
        dt=360
        K4=8.00*10**12
    else:
        print('Error: unsupported value of M. Only 42,63, 106, 170, and 213 are supported')
    
    #lmax=M
    # I = int(2*lmax + 1)#p.I 
    # J = int(lmax+1)#p.J

    
    lambdas=np.linspace(-np.pi, np.pi, num=I,endpoint=False) 
    [mus,w]=sp.roots_legendre(J)
    
    #mus_SH, w_SH = pysh.expand.SHGLQ(N)
    # Sine of the latitude
    #mus = mus_SH# p.mus #need negative to correspond to the transform
    # Weights for integrating
    #w = w_SH# p.w
        
    return N,I,J,dt,K4,lambdas,mus,w

def velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test):
    """Initializes the velocity components
    :param I: number of longitudes
    :type I: int
    :param J: number of latitudes
    :type J: int 

    :return: J by I data arrays for Uic and Vic
    :rtype: array of float64
    """
    Uic=np.full((J,I),0.0) #initialize
    Vic=np.full((J,I),0.0)
    
    if test==1:
        for i in range(I):
            for j in range(J):
                # Uic[j,i]=SU0*(np.cos(np.arcsin(mus[j]))*cosa +mus[j]*np.cos(lambdas[i])*sina)#assign according to init.f
                # Vic[j,i]=-SU0*np.sin(lambdas[i])*sina
                
                Uic[j,i]=SU0*(np.cos(np.arcsin(mus[j]))*cosa +mus[j]*np.cos(lambdas[i])*sina)*np.cos(np.arcsin(mus[j])) #assign according to init.f
                Vic[j,i]=-SU0*np.sin(lambdas[i])*sina*np.cos(np.arcsin(mus[j]))
    
    elif test==2: #Williamson Test 2
        for i in range(I):
            for j in range(J):   
                Uic[j,i] = SU0*(np.cos(np.arcsin(mus[j]))*cosa + np.cos(lambdas[i])*(mus[j])*sina)
                Vic[j,i] = -SU0*(np.sin(lambdas[i])*sina)
               

    elif test==10:
        for i in range(I):
            for j in range(J):
                Uic[j,i]=0#SU0*(np.cos(np.arcsin(mus[j]))*1 +mus[j]*np.cos(lambdas[i])*0)#assign according to init.f
                Vic[j,i]=0#-SU0*np.sin(lambdas[i])*0
    return Uic, Vic

def ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J):
    """Initializes the state auxiliary variables
    :param Uic: zonal velocity component
    :type Uic: array of complex128
    :param Vic: meridional velocity component
    :type Vic: array of complex128
    :param etaic0: initial eta
    :type etaic0:  array of complex128
    :param Phiic0: initial Phi
    :type Phiic0:  array of complex128
    :param mustile: reshaped mu array to fit the dimensions
    :type mustile:  array of complex128

    :return: J by I data arrays for eta0, eta1, delta0, delta1, and phi0, phi1
    :rtype:  array of complex128
    """
    Aic=np.multiply(Uic,etaic0) #A=U*\eta
    Bic=np.multiply(Vic,etaic0) #B=V*\eta
    Cic=np.multiply(Uic,Phiic0) #C=U*\Phi
    Dic=np.multiply(Vic,Phiic0) #D=V*\Phi
    
    Eic=np.zeros((J,I))
    #E=(U^2+V^2)/(2(1-mu^2))
    for i in range(I):
        for j in range(J):
            numerator=(Uic[j,i]**2+Vic[j,i]**2)
            denominator=(2*(1-mus[j]**2))
            Eic[j,i]=numerator/denominator

    return Aic, Bic, Cic, Dic, Eic

def coriolismn(M,omega):
    
    fmn=np.zeros([M+1,M+1]) 
    fmn[0,1]=omega/np.sqrt(0.375)
    
    return fmn