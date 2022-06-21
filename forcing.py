# -*- coding: utf-8 -*-
"""
This module contains the functions used for the evaluation of stellar forcing (insolation).
"""

import numpy as np 

def Phieqfun(Phibar,DPhieq,lambdas,mus,I,J,g):
    """
    Evaluates the equilibrium geopotential from Perez-Becker and Showman.

    Parameters
    ----------
    :param Phibar: 
        Mean geopotential
    :type Phibar: float
    
    :param DPhieq: 
        The difference between mean geopotential and the maximum geopotential
    :type DPhieq: float
    
    :param lambdas: 
        Uniformly spaced longitudes of length I.
    :type lambdas: array of float
    
    :param mus:
        Array of Gaussian latitudes of length J.
    :type mus: array of float

    :param I: number of longitudes    
    :type I: int
    
    :param J:
        number of latitudes.
    :type J: int

    :param g:
        Surface gravity, m/s^2.
    :type g: float

    Returns
    -------
    :return PhieqMat:  the equilibrium geopotential, array (J,I)
    :rtype: array of float

    """

    PhieqMat=Phibar*np.ones((J,I)) #initialize to flat nightside geopotentiAL
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (0,0)
            if  -np.pi/2<lambdas[i]<np.pi/2: #only force the dayside
                PhieqMat[j,i]=PhieqMat[j,i]+DPhieq*np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))     
    
    return PhieqMat



def Qfun(Phieq,Phi,Phibar,taurad):
    """
    Evaluates the radiative forcing on the geopotential. Corresponds to the 
    Q from Perez-Becker and Showman, but has an extra factor of g as we are 
    evaluating the geopotential, and they are evaluating the geopotential 
    height. 

    Parameters
    ----------
    :param Phieq: Equilibrium geopotential, (J,I)
    :type Phieq: array of float64
    
    :param Phi: Geopotential with the mean subtracted, (J,I)
    :type Phi: array of float64
    
    :param Phibar: Mean geopotential
    :type Phibar: float64
    
    :param taurad: radiative time scale, s
    :type taurad: float64


    Returns
    -------
    :return Q: Geopotential forcing, (J,I)
    :rtype: array of float64

    """
    #note Q is different from Perez-Becker and Showman by a factor of g (for consistency with Phi vs H)
    Q=(1/taurad)*(Phieq-(Phi+Phibar))

    return Q

def Qfun_with_rampup(Phieq,Phi,Phibar,taurad,t,dt):
    
    """
    Evaluates the radiative forcing on the geopotential, but slowly ramps up 
    the forcing to improve stability for short radiative timescales.

    Parameters
    ----------
    :param Phieq: Equilibrium geopotential, (J,I)
    :type Phieq: array of float64
    
    :param Phi: Geopotential with the mean subtracted, (J,I)
    :type Phi: array of float64
    
    :param Phibar: Mean geopotential
    :type Phibar: float64
    
    :param taurad: radiative time scale, s
    :type taurad: float64
    
    :param t: number of current timestep
    :type t: int
    
    :param dt: timestep length
    :type dt: float64


    Returns
    -------
    :return Q: Geopotential forcing, (J,I)
    :rtype: array of float64

    """
    
    #slowly ramp up
    if t*dt<15*3600:
        factor=t*dt/(15*3600)
    else:
        factor=1

    Q=factor*(1/taurad)*(Phieq-(Phi+Phibar))
    
    return Q


def Rfun(U,V,Q,Phi,Phibar, taudrag):
    """
    Evaluates the first and second component of the vector that expresses the 
    velocity forcing in Perez-Becker and Showman. The divergence and vorticity (F,G) 
    correspond to the forcing on the state variables delta and zeta, respectively.

    Parameters
    ----------
    :param U: zonal velocity component, (J,I)
    :type U: array of float
    
    :param V: meridional velocity component, (J,I)
    :type V: array of float
    
    :param Q: radiative forcing of geopotential, (J,I)
    :type Q: array of float
    
    :param Phi: geopotential with the mean subtracted, (J,I)
    :type Phi: array of float
    
    :param Phibar: mean geopotential
    :type Phibar: float
    
    :param taudrag: drag timescale,in seconds
    :type taudrag: float


    Returns
    -------
     :return:
        - F 
                Zonal component of the velocity forcing vector
        - G 
                Meridional component of the velocity forcing vector

    :rtype: array of float

    """
    
    Qclone=Q.copy()
    Qclone[Q<0]=0

    Ru=np.divide(np.multiply(-U,Qclone),Phi+Phibar)
    Rv=np.divide(np.multiply(-V,Qclone),Phi+Phibar)
    
    #reset to 0 if losing mass
    Ru[Q<0]=0
    Rv[Q<0]=0
    
     #if taudrag is infinity, only have the R componen   
    if taudrag!=-1:
        F=Ru-(U/taudrag)
        G=Rv-(V/taudrag)
        
    else:
        F=Ru
        G=Rv
        
    return F, G
