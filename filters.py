
"""
Created on Wed Jun 24 17:33:29 2020

@author: ek672
"""
import numpy as np

def modal_splitting(Xidataslice,alpha):
    """Applies the filter from Hack and Jacob (1992)
    :param Xidata: data array to be filtered
    :type Xidata: list
    :param alpha: filter coefficient
    :type alpha: float

    :return:
    :rtype: list
    """
    newxi=Xidataslice[1,:,:]+alpha*(Xidataslice[0]-2*Xidataslice[1]+Xidataslice[2])
    return newxi

def diffusion(Ximn,sigma):
    """ Applied the diffusion Filter described in Gelb and Gleeson (eq. 12)
    :param Ximn: the spectral coefficient data to be filtered
    :type Ximn: list
    :param currt: time t so that Ximndata(t) is filtered 
    :type currt: int 
    :param dt: time step
    :type t:
    :param a: radius of the planet, meters
    :type a: float
    :param K4: diffusion coefficient
    :type K4: float
    :param N:
    :type N: int
    """
    # sigma=np.zeros((M+1,N+1))
    # # nvec=np.arange(N+1)
    # # coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    # # sigmacoeff=(1+2*dt*K4*coeff/a**4)
    # # sigmas=np.divide(1,sigmacoeff)
    
    # nvec=np.arange(N+1)
    # coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    # const=K4*dt/a**4
    # sigmacoeff1=2*const*coeff
    # sigmacoeff=(1+sigmacoeff1)
    
    # #sigmacoeff=(1+2*p.dt*p.K4*coeff/p.a**4)
    # sigmas=np.divide(1,sigmacoeff)
    
    # for m in range(M+1):
    #     sigma[m,:]=sigmas
    
    newXimn=np.multiply(Ximn,sigma)
    return newXimn

# def diffusionPhi(Ximn,dt,a,K4,M, N):
#     """ Applied the diffusion Filter described in Gelb and Gleeson (eq. 12)
#     :param Ximn: the spectral coefficient data to be filtered
#     :type Ximn: list
#     :param currt: time t so that Ximndata(t) is filtered 
#     :type currt: int 
#     :param dt: time step
#     :type t:
#     :param a: radius of the planet, meters
#     :type a: float
#     :param K4: diffusion coefficient
#     :type K4: float
#     :param N:
#     :type N: int
#     """
#     sigma=np.zeros((M+1,N+1))
#     nvec=np.arange(N+1)
#     coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))
    
#     const=K4*dt/a**4
#     sigmacoeff1=2*const*coeff
#     sigmacoeff=(1+sigmacoeff1)
    
#     sigmas=np.divide(1,sigmacoeff)
#     for m in range(M+1):
#         sigma[m,:]=sigmas

#     newXimn=np.multiply(Ximn,sigma)
#     return newXimn

def sigma(M,N,K4,a,dt):
    sigma=np.zeros((M+1,N+1))
    # nvec=np.arange(N+1)
    # coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    # sigmacoeff=(1+2*dt*K4*coeff/a**4)
    # sigmas=np.divide(1,sigmacoeff)
    
    nvec=np.arange(N+1)
    coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    const=K4*dt/a**4
    sigmacoeff1=2*const*coeff
    sigmacoeff=(1+sigmacoeff1)
    
    #sigmacoeff=(1+2*p.dt*p.K4*coeff/p.a**4)
    sigmas=np.divide(1,sigmacoeff)
    
    for m in range(M+1):
        sigma[m,:]=sigmas
    
    return sigma

def sigmaPhi(M,N,K4,a,dt):
    sigma=np.zeros((M+1,N+1))
    nvec=np.arange(N+1)
    coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))
    
    const=K4*dt/a**4
    sigmacoeff1=2*const*coeff
    sigmacoeff=(1+sigmacoeff1)
    
    sigmas=np.divide(1,sigmacoeff)
    for m in range(M+1):
        sigma[m,:]=sigmas
    
    return sigma
    