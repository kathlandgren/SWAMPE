
"""
Created on Wed Jun  3 18:29:42 2020

@author: ek672

This file contains the initial conditions for the simulation"""

import numpy as np
import scipy.special as sp


def test1_init(a,omega,a1):
    """
    Initializes the parameters from Test 1 in Williamson et al. (1992),
    Advection of Cosine Bell over the Pole

    Parameters
    ----------
    :param a: 
        Planetary radius, in meters.
    :type a: float64
    
    :param omega: 
        Planetary rotation rate, in radians per second.
    :type omega: float64
    
    :param a1:
        Angle of advection, in radians.
    :type a1: float64

    Returns
    -------
    :return:
        - SU0 
            Amplitude parameter from Test 1 in Williamson et al. (1992)
        - sina 
            sine of the angle of advection.
        - cosa 
            cosine of the angle of advection.
        - etaamp 
            Amplitude of absolute vorticity.
        - Phiamp 
            Amplitude of geopotential.
    :rtype: float64

    """
    #Parameters for Test 1 in Williamson et al. (1992)
    SU0=2.0*np.pi*a/(3600.0*24*12) 
    sina=np.sin(a1) #sine of the angle of advection
    cosa=np.cos(a1)  #cosine of the angle of advection
    etaamp = 2.0*((SU0/a)+omega) #relative vorticity amplitude
    Phiamp = (SU0*a*omega + 0.5*SU0**2) #geopotential height amplitude
    
    return SU0, sina, cosa, etaamp, Phiamp

def state_var_init(I,J,mus,lambdas,test,etaamp,*args):
    """
    Initializes state variables.

    Parameters
    ----------
    :param I:
        number of longitudes.
    :type I: int
    :param J:
        number of latitudes.
    :type J: int
    :param mus:
        Array of Gaussian latitudes of length J.
    :type mus: array of float64
    :param lambdas: 
        Uniformly spaced longitudes of length I.
    :type lambdas: array of float64
    :param test:
        The number of the regime being tested.
    :type test: int
    :param etaamp:
        Amplitude of absolute vorticity.
    :type etaamp: float64
    *args : TYPE
        DESCRIPTION.

    Returns
    -------
    :return: 
        - etaic0 - Initial condition for absolute vorticity, (J,I).
        - etaic1 - Second initial condition for absolute vorticity, (J,I).
        - deltaic0 - Initial condition for divergence, (J,I).
        - deltaic1 - Second initial condition for divergence, (J,I).
        - Phiic0 - Initial condition for geopotential, (J,I).
        - Phiic1 - Second initial condition for geopotential, (J,I).

    :rtype: tuple of arrays of float64
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
    """
    Initializes the spectral parameters corresponding to the spectral resolution M

    Parameters
    ----------
    M : int
        Spectral truncation.

    Returns
    -------
    N : int
        DESCRIPTION.
    I : int
        number of longitudes.
    J : int
        number of latitudes.
    dt : float64
        length of time step.
    K4 : float64
        hyperviscosity parameter from Gelb and Gleeson.
    lambdas : array of float64
        Uniformly spaced longitudes of length I.
    mus : array of float64
        Array of Gaussian latitudes of length J.
    w : array of float 64
        Gaussian weights of length J.

    """
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
    """
    

    Parameters
    ----------
    :param I: number of latitudes.
    :type I:  int

    :param J: number of longitudes.
    :type J:  int
    
    :param SU0: Amplitude parameter from Test 1 in Williamson et al. (1992)
    :type SU0: float64

    :param cosa: cosine of the angle of advection.
    :type cosa: float64

    :param sina: sine of the angle of advection.
    :type sina: float64

    :param mus: Array of Gaussian latitudes of length J
    :type mus: array of float64
        
    :param lambdas: Array of uniformly spaces longitudes of length I.
    :type lambdas: array of float64
        
    :param test: The number of the regime being tested.
    :type test: int
        

    Returns
    -------
    :return: 
       - Uic - (J,I) array - the initial condition for the latitudinal velocity component,
       - Vic - (J,I) array - the initial condition for the meridional velocity component.
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
    """
    Initializes the Coriolis force in spectral space.
    
    Parameters
    ----------
    
    :param M: Spectral dimension.
    :type M: int
    
    :param omega: Planetary rotation rate, in radians per second.
    :type omega:  float64

    Returns
    -------
    :return: fmn The Coriolis force in spectral space of dimension (M+1, M+1).
    :rtype:  array of float64

    """
    
    fmn=np.zeros([M+1,M+1]) 
    fmn[0,1]=omega/np.sqrt(0.375)
    
    return fmn