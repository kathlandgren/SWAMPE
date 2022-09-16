# -*- coding: utf-8 -*-
"""
This module contains the functions associated with the modified Euler time-stepping scheme. 
The associated coefficients and the method are outlined in the Methods section of this documentation.
"""


import numpy as np
import scipy.special as sp
#local imports
from . import spectral_transform as st
from . import filters


## PHI tstep
def phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFm,dt,a,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
    """This function timesteps the geopotential Phi forward.

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
    :param test: number of test, defaults to None
    :type test: int
    :param t: number of current time step
    :type t: int

    :return: 
        - Phimntstep
                Updated spectral coefficients of geopotential
        - newPhitstep
                Updated geopotential
    :rtype: array of float
    """
    
    
    #use the "1" for state variables (the latest one)
    tstepcoeff1=tstepcoeff1/2
    tstepcoeff2=tstepcoeff2/2

    
    Phicomp1=st.fwd_leg(Phim1, J, M, N, Pmn, w)
    
    Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp2=st.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    Phicomp3prep=np.multiply(tstepcoeff1,Dm)
    Phicomp3=st.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
    Phicomp4=dt*Phibar*st.fwd_leg(deltam1, J, M, N, Pmn, w)
    
    
    deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp2=st.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
    

    deltacomp3prep=np.multiply(tstepcoeff1,Am)

    deltacomp3=st.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
    
    deltacomp5prep=np.multiply(tstepcoeff2,Em)
    deltacomp5=st.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    
    deltacomp5=np.multiply(narray,deltacomp5)
    
    Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))

    if forcflag==True:
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
    
        
        Phicomp1=st.fwd_leg(Phim1, J, M, N, Pmn, w)
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=st.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
        
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=st.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
        
        Phicomp4=dt*Phibar*st.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm+Fm)
        deltacomp2=st.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am-Gm)
    
        deltacomp3=st.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=st.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
       
        Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))

        Phiforcing=st.fwd_leg((dt)*PhiFm, J, M, N, Pmn, w)

        Phimntstep=Phimntstep+Phiforcing
    else:
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
    
        
        Phicomp1=st.fwd_leg(Phim1, J, M, N, Pmn, w)
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=st.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
        
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=st.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
        
        Phicomp4=dt*Phibar*st.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
        deltacomp2=st.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am)
    
        deltacomp3=st.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=st.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)

       
        Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))


    
    if diffflag==True:
        Phimntstep=filters.diffusion(Phimntstep, sigmaPhi) 
        
        
    
    newPhimtstep=st.invrs_leg(Phimntstep, I,J, M, N, Pmn)
    newPhitstep=st.invrs_fft(newPhimtstep, I)

    
    return Phimntstep,newPhitstep

def delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFm,dt,a,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
    """This function timesteps the divergence delta forward.
    
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
    :param test: number of test, defaults to None
    :type test: int
    :param t: number of current time step
    :type t: int

    :return: 
        - deltamntstep
                Updated spectral coefficients of divergence
        - newdeltatstep
                Updated divergence
    :rtype: array of float
    """
    
    tstepcoeff1=tstepcoeff1/2
    tstepcoeff2=tstepcoeff2/2
    
    
    deltacomp1=st.fwd_leg(deltam1, J, M, N, Pmn, w)
    
    deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp2=st.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
    

    deltacomp3prep=np.multiply(tstepcoeff1,Am)

    deltacomp3=st.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
    
    deltacomp4prep=np.multiply(tstepcoeff2,Phim1)
    deltacomp4=st.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
    
    deltacomp4=np.multiply(narray,deltacomp4)
    
    deltacomp5prep=np.multiply(tstepcoeff2,Em)
    deltacomp5=st.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    
    deltacomp5=np.multiply(narray,deltacomp5)
    
    Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp2=st.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)


    Phicomp3prep=np.multiply(tstepcoeff1,Dm)
    Phicomp3=st.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)

    deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2


    if forcflag==True:
        
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
        
        
        deltacomp1=st.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm+Fm)
        deltacomp2=st.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am-Gm)
    
        deltacomp3=st.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        deltacomp4prep=np.multiply(tstepcoeff2,Phim1)
        deltacomp4=st.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
        
        deltacomp4=np.multiply(narray,deltacomp4)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=st.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
        
    
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=st.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=st.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
        deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2

        Phiforcing=np.multiply(narray,st.fwd_leg((dt/2)*PhiFm, J, M, N, Pmn, w))/a**2

        deltamntstep=deltamntstep+Phiforcing
        
    else:
            
        
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
        
        
        deltacomp1=st.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm+Fm)
        deltacomp2=st.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am-Gm)
    
        deltacomp3=st.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        deltacomp4prep=np.multiply(tstepcoeff2,Phim1)
        deltacomp4=st.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
        
        deltacomp4=np.multiply(narray,deltacomp4)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=st.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
        
    
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=st.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=st.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
     
        deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2
   
 
        
        

    
    if diffflag==True:
        deltamntstep=filters.diffusion(deltamntstep, sigma)

    newdeltamtstep=st.invrs_leg(deltamntstep, I,J, M, N, Pmn)
    newdeltatstep=st.invrs_fft(newdeltamtstep, I)
    return deltamntstep,newdeltatstep


def eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFm,dt,a,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
     """This function timesteps the absolute vorticity eta forward.

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
    :param test: number of test, defaults to None
    :type test: int
    :param t: number of current time step
    :type t: int

    :return: 
        - etamntstep
                Updated spectral coefficients of absolute vorticity
        - newetatstep
                Updated absolute vorticity
    :rtype: array of float
    """
    
    if forcflag==True:
        etacomp1=st.fwd_leg(etam1, J, M, N, Pmn, w)
    
        etacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Am-Gm)
        etacomp2=st.fwd_leg(etacomp2prep, J, M, N, Pmn, w)
        
        etacomp3prep=np.multiply(tstepcoeff1,Bm+Fm)
        etacomp3=st.fwd_leg(etacomp3prep, J, M, N, Hmn, w)
        etamntstep=etacomp1-etacomp2+etacomp3
        
    else:
        tstepcoeff1=tstepcoeff1/2

        etacomp1=st.fwd_leg(etam1, J, M, N, Pmn, w)
    
        etacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Am)
        etacomp2=st.fwd_leg(etacomp2prep, J, M, N, Pmn, w)
        
        etacomp3prep=np.multiply(tstepcoeff1,Bm)
        etacomp3=st.fwd_leg(etacomp3prep, J, M, N, Hmn, w)
    
    
    etamntstep=etacomp1-etacomp2+etacomp3

    
    if diffflag==True:
        etamntstep=filters.diffusion(etamntstep, sigma)
    
    newetamtstep=st.invrs_leg(etamntstep, I,J, M, N, Pmn)
    newetatstep=st.invrs_fft(newetamtstep, I)
    return etamntstep,newetatstep
