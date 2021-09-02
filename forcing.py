# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:28:45 2020

@author: ek672

This module contains the functions used for the evaluation of forcing.
"""

import numpy as np 

import testing_plots
import scipy.special as sp

def Phieqfun(Phibar,DPhieq,lambdas,mus,I,J,g):
    """
    Evaluates the equilibrium geopotential from Perez-Becker and Showman.

    Parameters
    ----------
    :param Phibar: 
        Mean geopotential
    :type Phibar: float64
    
    :param DPhieq: 
        The difference between mean geopotential and the maximum geopotential
    :type DPhieq: float64
    
    :param lambdas: 
        Uniformly spaced longitudes of length I.
    :type lambdas: array of float64 
    
    :param mus:
        Array of Gaussian latitudes of length J.
    :type mus: array of float64

    :param I: number of longitudes    
    :type I: int
    
    :param J:
        number of latitudes.
    :type J: int

    :param g:
        Surface gravity, m/s^2.
    :type g: float64

    Returns
    -------
    :return PhieqMat:  the equilibrium geopotential, array (J,I)
    :rtype: array of float64

    """

    PhieqMat=Phibar*np.ones((J,I)) #initialize to flat nightside temperature
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (0,0)
            if  -np.pi/2<lambdas[i]<np.pi/2:
                PhieqMat[j,i]=PhieqMat[j,i]+DPhieq*np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))     

    return PhieqMat


def insolation(L_star,distance,sigmaSB):
    insolation=((L_star)/(4*np.pi*(distance**2)))/sigmaSB
    
    return insolation

# def Qfun(heq,Phi,Phibar,taurad,g):
#     Q=(1/taurad)*(heq-(Phi+Phibar)/g)
#     #Q=(1/taurad)*(heq-(Phi)/g)
#     return Q
def DoubleGrayTEqfun(Phibar,insolation,lambdas,mus,I,J,k1,k2,p,g,R,Cp,sigma):
    """
    Evaluates the equilibrium temperature from Langton and Laughlin (2008).
    
    Parameters
    ----------
    :param Phibar: 
        Mean geopotential
    :type Phibar: float64
    :param DPhieq: 
        The difference between mean geopotential and the maximum geopotential
    :type DPhieq: float64
    
    :param lambdas: 
        Uniformly spaced longitudes of length I.
    :type lambdas: array of float64 
    
    :param mus:
        Array of Gaussian latitudes of length J.
    :type mus: array of float64
    
    :param I: number of longitudes    
    :type I: int
    
    :param J:
        number of latitudes.
    :type J: int
    
    :param k1: visible opacity
    :type k1: float64
    
    :param k2: infra-red opacity
    :type k2: float64
    
    :param p: pressure
    :type p: float64
    
    :param g:
        Surface gravity, m/s^2.
    :type g: float64
    
    :param R: specific gas constant
    :type R: float64
        
    :param Cp: heat capacity at constant pressure
    :type Cp: float64
    
    :param sigma: Stefan-Boltmann constant
    :type sigma: float64
    
    Returns
    -------
    :return TeqMat:  the equilibrium temperature, array (J,I)
    :rtype: array of float64
    
    """

    TeqMat=((Phibar/R)**4)*np.ones((J,I))

    x=np.exp(-k1*p/g)
    
    #insolation=DPhieq#(DPhieq/R)**4
    ratio=(k1/k2)
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (0,0)
            if  -np.pi/2<lambdas[i]<np.pi/2:
                
                
                ## implementation of Langton formulation
                # ss_angle_sec=(1/(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))))
                
                # day_forcing=ratio*(insolation*(x**ss_angle_sec))
                
                # TeqMat[j,i]=TeqMat[j,i]+day_forcing
                
                
                ## cosine bell experiment
                
                ss_angle=(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2)))
                
                day_forcing=ratio*insolation*(x*ss_angle)
                
                TeqMat[j,i]=TeqMat[j,i]+day_forcing
                
                #TeqMat[j,i]=TeqMat[j,i]+(k1/k2)*(((DPhieq/R)**4)*(x**(1/(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))))))#+(Phibar/R)**4
                #PhieqMat[j,i]=PhieqMat[j,i]+k1*DPhieq*(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2)))*x**(1/(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))))        

    return TeqMat

def DoubleGrayPhiForcing(TeqMat,Phidata,Phibar,k2,sigma,Cp,R):
    
    """
    Evaluates the radiative forcing on the geopotential using the double-gray 
    scheme. Corresponds to the radiative forcing from Langton and Laughlin
    (2008), converted to geopotential from temperature using the ideal gas law
    relationship.


    Parameters
    ----------
    :param TeqMat: Equilibrium temperature, (J,I)
    :type TeqMat: array of float64
    
    :param Phidata: Geopotential with the mean subtracted, (J,I)
    :type Phidata: array of float64
    
    :param Phibar: Mean geopotential
    :type Phibar: float64
    
    :param k2: infra-red opacity
    :type k2: float64

    :param sigma: Stefan-Boltzmann constant, W/(m^2 K^4)
    :type sigma: float64

    Returns
    -------
    :return Q: Geopotential forcing, (J,I)
    :rtype: array of float64

    """
    outer_coeff=(sigma*k2*R/Cp)
    #sci_comp_step=(TeqMat-((Phidata+Phibar)/R)**4)/10**11
    A=(TeqMat-((Phidata+Phibar)/R)**4)
    # print(np.min(TeqMat))
    # print(np.min(((Phidata+Phibar)/R)**4))
    Q=outer_coeff*A
    # print('Max Q is '+str(np.max(Q)))
    
    
    
    # I=192
    # J=96
    
    # lambdas=np.linspace(-np.pi, np.pi, num=I,endpoint=False) 
    # [mus,w]=sp.roots_legendre(J)
    # testing_plots.physical_plot(TeqMat,mus,lambdas)
    # testing_plots.physical_plot(((Phidata+Phibar)/R)**4,mus,lambdas)
    # testing_plots.physical_plot(Q,mus,lambdas)
    return Q
    
    

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
    #note Q is different from Perez-Becker and Showman, our Q is PBS-Q*g
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
    #note Q is different from Perez-Becker and Showman, our Q is PBS-Q*g
    
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
    :type U: array of float64
    
    :param V: meridional velocity component, (J,I)
    :type V: array of float64
    
    :param Q: radiative forcing of geopotential, (J,I)
    :type Q: array of float64
    
    :param Phi: geopotential with the mean subtracted, (J,I)
    :type Phi: array of float64
    
    :param Phibar: mean geopotential
    :type Phibar: float64
    
    :param taudrag: drag timescale,in seconds
    :type taudrag: float64


    Returns
    -------
     :return:
        - F 
                First component of the velocity forcing vector
        - G 
                Second component of the velocity forcing vector

    :rtype: array of float64

    """
    
    Qtest=Q.copy()
    Qtest[Q<0]=0
    # Ru=np.divide(np.multiply(-U,Q),Phi+Phibar)
    # Rv=np.divide(np.multiply(-V,Q),Phi+Phibar)
    Ru=np.divide(np.multiply(-U,Qtest),Phi+Phibar)
    Rv=np.divide(np.multiply(-V,Qtest),Phi+Phibar)
    
    #reset to 0 if losing mass
    Ru[Q<0]=0
    Rv[Q<0]=0
    
     #if taudrag is infinity, only have the R component from Perez-Becker and Showman    
    if taudrag!=-1:
        F=Ru-(U/taudrag)
        G=Rv-(V/taudrag)
        
    else:
        F=Ru
        G=Rv
        
    return F, G
