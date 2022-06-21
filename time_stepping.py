# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:03:58 2020

@author: ek672
"""

import numpy as np

#local imports 
import params as p
import spectral_transform as st


expflag=p.expflag
if expflag==1:
    import explicit_tdiff as tdiff
else:

    import modEuler_tdiff as tdiff


def tstepping(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,marray,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
    """
    
    :param etam0: DESCRIPTION
    :type etam0: TYPE
    :param etam1: DESCRIPTION
    :type etam1: TYPE
    :param deltam0: DESCRIPTION
    :type deltam0: TYPE
    :param deltam1: DESCRIPTION
    :type deltam1: TYPE
    :param Phim0: DESCRIPTION
    :type Phim0: TYPE
    :param Phim1: DESCRIPTION
    :type Phim1: TYPE
    :param I: DESCRIPTION
    :type I: TYPE
    :param J: DESCRIPTION
    :type J: TYPE
    :param M: DESCRIPTION
    :type M: TYPE
    :param N: DESCRIPTION
    :type N: TYPE
    :param Am: DESCRIPTION
    :type Am: TYPE
    :param Bm: DESCRIPTION
    :type Bm: TYPE
    :param Cm: DESCRIPTION
    :type Cm: TYPE
    :param Dm: DESCRIPTION
    :type Dm: TYPE
    :param Em: DESCRIPTION
    :type Em: TYPE
    :param Fm: DESCRIPTION
    :type Fm: TYPE
    :param Gm: DESCRIPTION
    :type Gm: TYPE
    :param Um: DESCRIPTION
    :type Um: TYPE
    :param Vm: DESCRIPTION
    :type Vm: TYPE
    :param fmn: DESCRIPTION
    :type fmn: TYPE
    :param Pmn: DESCRIPTION
    :type Pmn: TYPE
    :param Hmn: DESCRIPTION
    :type Hmn: TYPE
    :param w: DESCRIPTION
    :type w: TYPE
    :param tstepcoeff: DESCRIPTION
    :type tstepcoeff: TYPE
    :param tstepcoeff2: DESCRIPTION
    :type tstepcoeff2: TYPE
    :param tstepcoeffmn: DESCRIPTION
    :type tstepcoeffmn: TYPE
    :param marray: DESCRIPTION
    :type marray: TYPE
    :param mJarray: DESCRIPTION
    :type mJarray: TYPE
    :param narray: DESCRIPTION
    :type narray: TYPE
    :param PhiFM: DESCRIPTION
    :type PhiFM: TYPE
    :param dt: DESCRIPTION
    :type dt: TYPE
    :param a: DESCRIPTION
    :type a: TYPE
    :param K4: DESCRIPTION
    :type K4: TYPE
    :param Phibar: DESCRIPTION
    :type Phibar: TYPE
    :param taurad: DESCRIPTION
    :type taurad: TYPE
    :param taudrag: DESCRIPTION
    :type taudrag: TYPE
    :param forcflag: DESCRIPTION
    :type forcflag: TYPE
    :param diffflag: DESCRIPTION
    :type diffflag: TYPE
    :param sigma: DESCRIPTION
    :type sigma: TYPE
    :param sigmaPhi: DESCRIPTION
    :type sigmaPhi: TYPE
    :param test: DESCRIPTION
    :type test: TYPE
    :param t: DESCRIPTION
    :type t: TYPE
    
    
    :return: DESCRIPTION
    :rtype: TYPE

    """
    
    newPhimn,newPhitstep2=tdiff.phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t)
    newdeltamn,newdeltatstep2=tdiff.delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t)
    newetamn,newetatstep2=tdiff.eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t)
    
    Unew,Vnew=st.invrsUV(newdeltamn,newetamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    
    return newetamn,newetatstep2,newdeltamn,newdeltatstep2,newPhimn,newPhitstep2,Unew,Vnew


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
    :type dt: float64
    
    :param a: planetary radius, m
    :type a: float64

    :return tstepcoeff2: time stepping coefficients of size (J,M+1)
    :rtype: array of float64

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
    :rtype: array of float64

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
    :type mus: array of float64
    
    :param a: planetary radius, m
    :type a: float64
    
    :return tstepcoeff: coefficients 2dt/(a(1-mus^2)) in a matrix of size (J,M+1)
    :rtype: array of float64
    
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
    :rtype: array of float64
    """
    
    mJarray=np.zeros((J,M+1))
    mtemp=np.arange(M+1)
    for j in range(J):
        mJarray[j,:]=mtemp
    return mJarray

def marray(M,N):
    """
    Computes coefficients equal to m=0,1,...,M

   
    :param M: spectral dimension
    :type M: int   

    :param N: highest degreee of associated Legendre polynomials
    :type N: int
    
    :return marray: coefficient m in a matrix of size (M+1, N+1)
    :rtype: array of float64
   
    """
    marray=np.zeros((M+1,N+1)) #TODO make this an input, compute once
    mtemp=np.arange(M+1)
    for n in range(N+1):
        marray[:,n]=mtemp
    return marray