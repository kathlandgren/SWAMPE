# -*- coding: utf-8 -*-
"""
This module contains the function that calls an explicit or an implicit time-stepping scheme,
as well as the functions that compute the arrays of coefficients involved in time-stepping. 
"""

import numpy as np

#local import
from . import spectral_transform as st




def tstepping(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,marray,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,expflag,sigma,sigmaPhi,test,t):
    """
    
    Calls the timestepping scheme.
    
    :param etam0: Fourier coefficents of absolute vorticity for one time step
    :type etam0: array of float
    :param etam1: Fourier coefficents of absolute vorticity for the following time step
    :type etam1: array of float
    :param deltam0: Fourier coefficents of divergence for one time step
    :type deltam0: array of float
    :param deltam1: Fourier coefficents of divergence for the following time step
    :type deltam1:  array of float
    :param Phim0: Fourier coefficents of geopotential for one time step
    :type Phim0: array of float
    :param Phim1: Fourier coefficents of geopotential for the following time step
    :type Phim1: array of float
    :param I: number of longitudes
    :type I: int
    :param J: number of Gaussian latitudes
    :type J: int
    :type M: int
    :param N: highest degree of the Legendre functions for m=0
    :type N: int
    :param Am: Fourier coefficients of the nonlinear component A=U*eta
    :type Am: array of float
    :param Bm: Fourier coefficients of the nonlinear component B=V*eta
    :type Bm: array of float
    :param Cm: Fourier coefficients of the nonlinear component C=U*Phi
    :type Cm: array of float
    :param Dm: Fourier coefficients of the nonlinear component D=V*Phi
    :type Dm: array of float
    :param Em: Fourier coefficients of the nonlinear component E=(U^2+V^2)/(2(1-mu^2))
    :type Em: array of float
    :param Fm: Fourier coefficients of the zonal component of wind forcing
    :type Fm: array of float
    :param Gm:  Fourier coefficients of the meridional component of wind forcing
    :type Gm: array of float
    :param Um: Fourier coefficients of the zonal component of wind
    :type Um: array of float
    :param Vm: Fourier coefficients of the meridional component of wind
    :type Vm: array of float
    :param fmn: spectral coefficients of the Coriolis force
    :type fmn: array of float
    :param Pmn: associated legendre functions evaluated at the Gaussian latitudes mus  up to wavenumber M
    :type Pmn: array of float
    :param Hmn: derivatives of the associated legendre functions evaluated at the Gaussian latitudes mus  up to wavenumber M
    :type Hmn: array of float
    :param w: Gauss Legendre weights
    :type w: array of float
    
    :param tstepcoeff: coefficient for time-stepping of the form 2dt/(a(1-mus^2))
    :type tstepcoeff: array of float
    
    :param tstepcoeff2: time stepping coefficient of the form 2dt/a^2
    :type tstepcoeff2: array of float
    
    :param tstepcoeffmn: an array of coefficients a/(n(n+1))
    :type tstepcoeffmn: array of float
    
    :param marray: coefficients equal to m=0,1,...,M in a matrix M+1xN+1
    :type marray: array of float
    
    :param mJarray:  coefficients equal to m=0,1,...,M in a matrix M+1xJ
    :type mJarray: array of float
    
    :param narray: array n(n+1) in a matrix M+1xN+1
    :type narray: array of float
    :param PhiFm: Fourier coefficients of the geopotential forcing 
    :type PhiFM: array of float
    :param dt: time step, in seconds
    :type dt: float
    :param a: planetary radius, m
    :type a: float
    :param K4: fourth-degree hyperdiffusion filter coefficient
    :type K4: float
    :param Phibar: time-invariant spatial mean geopotential, height of the top layer
    :type Phibar: float
    :param taurad: radiative timescale
    :type taurad: float
    :param taudrag: drag timescale
    :type taudrag: float
    :param forcflag: forcing flag
    :type forcflag: float
    
    :param diffflag: hyperdiffusion filter flag
    :type diffflag: float
    
    :param sigma: hyperdiffusion filter coefficients for absolute vorticity and divergence
    :type sigma: array of float
    :param sigmaPhi: hyperdiffusion filter coefficients for geopotential
    :type sigmaPhi: array of float
    :param test: number of test
    :type test: TYPE
    :param t: number of current time step
    :type t: int
    
    
    :return: 
        - newetamn 
                Updated spectral coefficients of absolute vorticity
        - newetatstep
                Updated absolute vorticity
        - newdeltamn
                Updated spectral coefficients of divergence
        - newdeltatstep
                Updated divergence
        - newPhimn
                Updated spectral coefficients of geopotential
        - newPhitstep
                Updated geopotential
        - Unew
                Updated zonal winds
        - Vnew 
                Updated meridional winds
                
    :rtype: array of float

    """
    
    #import explicit or implicit time difference scheme
    if expflag==1:
        from . import explicit_tdiff as tdiff
    else:
        from . import modEuler_tdiff as tdiff

    
    newPhimn,newPhitstep=tdiff.phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t)
    newdeltamn,newdeltatstep=tdiff.delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t)
    newetamn,newetatstep=tdiff.eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t)
    
    Unew,Vnew=st.invrsUV(newdeltamn,newetamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
    return newetamn,newetatstep,newdeltamn,newdeltatstep,newPhimn,newPhitstep,Unew,Vnew


def tstepcoeffmn(M,N,a):
    """
    Generates the coefficient that multiplies spectral components in invrsUV
    :param M: highest wave number
    :type M: int
    :param N: highest degree of associated legendre polynomial for m=0
    :param a: radius of the planet, in meters
    :type a: float
    :param ns: M+1xN+1xIxJ array that only varies along the N axis
    :type ns: list
    
    :return: an array of coefficients a/(n(n+1))
    :rtype: list
    """
    coeff=np.multiply(np.arange(N+1),np.arange(N+1)+1)
    coeff[0]=1.
    tstepcoeff=a/coeff
    tstepcoeff[0]=0
    tstepcoeffmn=np.zeros([M+1,N+1])

    for m in range(M+1):
        tstepcoeffmn[m,:]=tstepcoeff
        
    return tstepcoeffmn


def tstepcoeff2(J,M,dt,a):
    """
    Computes the time stepping coefficient of the form 2dt/a^2 from Hack and 
    Jakob (1992)
    
    :param J: number of Gaussian latitudes
    :type J: int
    
    :param M: spectral dimension
    :type M: int
    
    :param dt: time step length, s
    :type dt: float
    
    :param a: planetary radius, m
    :type a: float

    :return tstepcoeff2: time stepping coefficients of size (J,M+1)
    :rtype: array of float

    """
    tstepcoeff2=np.zeros((J,M+1))
    for m in range(M+1):
        tstepcoeff2[:,m]=(2*dt/a**2)
    return tstepcoeff2
    

def narray(M,N):
    """
    Computes the array n(n+1).

    :param M: spectral dimension
    :type M: int
    
    :param N: highest degreee of associated Legendre polynomials
    :type N: int
    
    :return narray: coefficients n(n+1) in a matrix of size (M+1,N+1)
    :rtype: array of float

    """
    narray=np.zeros((M+1,N+1))
    for m in range(M+1):
        for n in range(N+1):
            narray[m,n]=n*(n+1)
    return narray

def tstepcoeff(J,M,dt,mus,a):
    """
    Computes a coefficient for time-stepping of the form 2dt/(a(1-mus^2))
    from Hack and Jakob (1992)
    
    :param J: number of Gaussian latitudes
    :type J: int
    
    :param M: spectral dimension
    :type M: int
    
    :param mus: array of Gaussian latitudes of length J
    :type mus: array of float
    
    :param a: planetary radius, m
    :type a: float
    
    :return tstepcoeff: coefficients 2dt/(a(1-mus^2)) in a matrix of size (J,M+1)
    :rtype: array of float
    
    """
    tstepcoeff=np.zeros((J,M+1))
    for m in range(M+1):
       tstepcoeff[:,m]=(2*dt)/(a*(1-mus**2))
    return tstepcoeff

def mJarray(J,M):
    """
    Computes coefficients equal to m=0,1,...,M

    :param J: number of Gaussian latitudes
    :type J: int
    
    :param M: spectral dimension
    :type M: int   

    :return mJarray: coefficient m in a matrix of size (J, M+1)
    :rtype: array of float
    """
    
    mJarray=np.zeros((J,M+1))
    mtemp=np.arange(M+1)
    for j in range(J):
        mJarray[j,:]=mtemp
    return mJarray

def marray(M,N):
    """
    Computes coefficients equal to m=0,1,...,M

   
    :param M: highest wavenumber of associated Legendre polynomials
    :type M: int   

    :param N: highest degreee of associated Legendre polynomials
    :type N: int
    
    :return marray: coefficient m in a matrix of size (M+1, N+1)
    :rtype: array of float
   
    """
    marray=np.zeros((M+1,N+1)) #TODO make this an input, compute once
    mtemp=np.arange(M+1)
    for n in range(N+1):
        marray[:,n]=mtemp
    return marray