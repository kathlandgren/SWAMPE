
"""
Created on Wed Jun  3 18:29:42 2020

@author: ek672

This file contains the initial conditions for the simulation
"""

import numpy as np
import scipy.special as sp

#local imports
import params as p
import pyshtools as pysh


##specifies the testing regime: 
# 1 -- test 1 from Williamson (advection of cosine bell)
# 2 -- test 2 from Williamson (global steady state nonlinear zonal geostrophic flow with compact support)
# 3 -- test 3 from Williamson (advection of cosine bell)
# 6 -- test 6 from Williamson (Rossby-Haurwitz wave)
# 10 -- Hot Jupiter (PBS) 


def test1_init(a,omega,a1):
    #Parameters for Test 1 in Williamson et al. (1992)
    SU0=2.0*np.pi*a/(3600.0*24*12) 
    sina=np.sin(a1) #sine of the angle of advection
    cosa=np.cos(a1)  #cosine of the angle of advection
    etaamp = 2.0*((SU0/a)+omega) #relative vorticity amplitude
    Phiamp = (SU0*a*omega + 0.5*SU0**2) #geopotential height amplitude
    
    return SU0, sina, cosa, etaamp, Phiamp
    

def state_var_init(I,J,mus,lambdas,g,omega,Phibar,test,*args):
    """Initializes the state variables
    :param I: number of longitudes
    :type I: int
    :param J: number of latitudes
    :type J: int 

    :return: J by I data arrays for eta0, eta1, delta0, delta1, and phi0, phi1
    :rtype: arrays of float64
    """

    etaic0=np.zeros((J,I))    
    zetaic0=np.zeros((J,I))
    Phiic0=np.zeros((J,I))
    deltaic0=np.zeros((J,I))
    
    if test<3:
        a,sina,cosa,etaamp,Phiamp,f_latlon=args
    
    if test==1: # Williamson Test 1 as documented in stswm FORTRAN implementation (see init.i)
        bumpr=a/3 #radius of the bump
        
        #starting coordinates for the cosine bell to be advected
        mucenter=0
        lambdacenter=3*np.pi/2
        for i in range(I):
            for j in range(J):
                etaic0[j,i]=etaamp*(-np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*sina+(mus[j])*cosa)
               
                dist=p.a*np.arccos(mucenter*mus[j]+np.cos(np.arcsin(mucenter))*np.cos(np.arcsin(mus[j]))*np.cos(lambdas[i]-lambdacenter))
                if dist < bumpr:
                    Phiic0[j,i]=(Phibar/2)*(1+np.cos(np.pi*dist/(bumpr)))
            
        #zeta since swarztrauber timesteps zeta, and not eta
        zetaic0=etaic0-f_latlon
                
    
    if test==2: #Williamson Test 2 as documented in stswm FORTRAN implementation (see init.i)
        for i in range(I):
            for j in range(J):
                latlonarg = -np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*sina+(mus[j])*cosa
                etaic0[j,i]=etaamp*(latlonarg)

                Phiic0[j,i]=((Phibar-Phiamp)*(latlonarg)**2)/g
                
        #zeta since Swarztrauber timesteps zeta, and not eta
        zetaic0=etaic0-f_latlon
    
    elif test==10: #PBS Hot Jupiter
        for i in range(I):
            for j in range(J):
                
                zetaic0[j,i]=0#etaamp*(-np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*0+(mus[j])*1)
                
               
    zetaic1=zetaic0 #need two time steps to initialize

    deltaic1=deltaic0

    Phiic0=Phiic0+Phibar
    Phiic1=Phiic0
    
    return zetaic0, zetaic1, deltaic0, deltaic1, Phiic0, Phiic1

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
    

    I=I-1 #reset for SHtools calculations
    #set longitude array
    lambdas=np.linspace(0, 2*np.pi, num=I,endpoint=False) 
    
    #set sin(latitude) array and the corresponding Gaussian weights
    mus, w = pysh.expand.SHGLQ(J-1)
    
    #normalization for the spherical harmonics
    normnum = 1
    
    return N,I,J,dt,K4,lambdas,mus,w,normnum

def velocity_init(I,J,mus,lambdas,test, *args):
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
    
    if test<3:
        SU0,cosa,sina=args
    
    if test==1: #Williamson Test 1
        for i in range(I):
            for j in range(J):
                Uic[j,i]=SU0*(np.cos(np.arcsin(mus[j]))*cosa +mus[j]*np.cos(lambdas[i])*sina)#assign according to init.f
                Vic[j,i]=-SU0*np.sin(lambdas[i])*sina
                
                #Uic[j,i]=SU0*(np.cos(np.arcsin(mus[j]))*cosa +mus[j]*np.cos(lambdas[i])*sina)*np.cos(np.arcsin(mus[j])) #assign according to init.f
                #Vic[j,i]=-SU0*np.sin(lambdas[i])*sina*np.cos(np.arcsin(mus[j]))
    
    elif test==2: #Williamson Test 2
        for i in range(I):
            for j in range(J):   
                Uic[j,i] = SU0*(np.cos(np.arcsin(mus[j]))*cosa + np.cos(lambdas[i])*(mus[j])*sina)
                Vic[j,i] = -SU0*(np.sin(lambdas[i])*sina)
   
    
    elif test==10: #PBS Hot Jupiter
        for i in range(I):
            for j in range(J):
                Uic[j,i]=0#SU0*(np.cos(np.arcsin(mus[j]))*cosa +mus[j]*np.cos(lambdas[i])*sina)#assign according to init.f
                Vic[j,i]=0#-SU0*np.sin(lambdas[i])*sina
    return Uic, Vic


def f_latlon(mus,lambdas,I,J,omega,a1,test):
    fMAT=np.zeros((J,I))
    for i in range(I):
        if test==1: #Williamson Test 1
            fMAT[:,i]=mus
        elif test==2: #Williamson Test 2, so that the flow can be specified with the spherical coordinate poles no necessarily coincident with the roation axis
            fMAT[:,i]=mus      
        # for j in range(J):
            #     fMAT[j,i]=-np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*np.sin(a1)+(mus[j])*np.cos(a1) #eq. 96 in Williamson et al. (1992) p.218
        elif test==10:
            fMAT[:,i]=mus            
    fMAT=fMAT*2*omega
    return fMAT