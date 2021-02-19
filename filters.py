
"""
Created on Wed Jun 24 17:33:29 2020

@author: ek672
"""
import numpy as np
import mpmath
import pyshtools as pysh

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
    print(np.mean(Xidataslice[1:,:]))
    print(np.mean(alpha*(Xidataslice[0]-2*Xidataslice[1]+Xidataslice[2])))
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
    
    newXimn=np.multiply(Ximn,sigma)
    return newXimn


def sigma(M,N,K4,a,dt):
    sigma=np.zeros((M+1,N+1))
    
    nvec=np.arange(N+1)
    coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    const=K4*dt/a**4
    sigmacoeff1=2*const*coeff
    sigmacoeff=(1+sigmacoeff1)
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


def hyperviscfun(statevar,w,mus,J,M,musMAT,sigma,normnum,lmax):
    """Computes the hyperviscosity given sigma and state variable
            
        :param M: 
        :type M: 
        :param mus: 
        :type mus:
        
        :return: 
        :rtype: 
        :return: 
        :rtype: 
        """
        
    #lmax=M+1, lmax-1=M
    svlm = pysh.expand.SHExpandGLQ(statevar, w, mus, norm=normnum,csphase=1,lmax_calc=lmax-1)
    
    svlm[0,:,:]=np.multiply(np.transpose(sigma),svlm[0,:,:])
    svlm[1,:,:]=np.multiply(np.transpose(sigma),svlm[1,:,:])
    
    svlmPad=np.zeros((2,J,J))
    svlmPad[0,:M+1,:M+1]=svlm[0,:,:]
    svlmPad[1,:M+1,:M+1]=svlm[1,:,:]
    
    hypervisc1=pysh.expand.MakeGridGLQ(svlmPad, mus, norm=normnum)

    return hypervisc1


##compute lat-lon hyperviscosity
# def hyperviscfun(statevar,w,mus,J,M,musMAT,nMAT2,normnum,lmax):
#     """Computes the Laplacian according to step 7.5 in Swarztrauber (1996)
            
#         :param M: 
#         :type M: 
#         :param mus: 
#         :type mus:
        
#         :return: 
#         :rtype: 
#         :return: 
#         :rtype: 
#         """
        
# ##Compute two Laplacians
#     #lmax=M+1, lmax-1=M
#     svlm = pysh.expand.SHExpandGLQ(statevar, w, mus, norm=normnum,csphase=1,lmax_calc=lmax-1)
    
#     nsqMAT=np.multiply(nMAT2,nMAT2)
    
#     svlm[0,:,:]=np.multiply(np.transpose(nsqMAT),svlm[0,:,:])
#     svlm[1,:,:]=np.multiply(np.transpose(nsqMAT),svlm[1,:,:])
    
#     svlmPad=np.zeros((2,J,J))
#     svlmPad[0,:M+1,:M+1]=svlm[0,:,:]
#     svlmPad[1,:M+1,:M+1]=svlm[1,:,:]
    
#     hypervisc1=pysh.expand.MakeGridGLQ(svlmPad, mus, norm=normnum)

#     return hypervisc1

def hyperviscfilter(delta,zeta,Phi,a,K4,w,mus,J,M,musMAT,sigma,sigmaPhi,normnum,lmax):
    
    deltaHV=hyperviscfun(delta,w,mus,J,M,musMAT,sigma,normnum,lmax)#-4*delta
    zetaHV=hyperviscfun(zeta,w,mus,J,M,musMAT,sigma,normnum,lmax)#-4*zeta #zeta might need to be zeta+f_latlon
    PhiHV=hyperviscfun(delta,w,mus,J,M,musMAT,sigmaPhi,normnum,lmax)    
    
    # deltaHV=-deltaHV1*K4#/a**4
    # zetaHV=-zetaHV1*K4#/a**4
    # PhiHV=-PhiHV1*K4 # /a**4   
    
    return deltaHV,zetaHV,PhiHV
    
    
    
    
    
## Spectral Viscosity
def sigmaSV(M,N,qvec,a,dt):
    sigma=np.array(np.zeros((M+1,N+1)),dtype=np.longfloat)
    # nvec=np.arange(N+1)
    # coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    # sigmacoeff=(1+2*dt*K4*coeff/a**4)
    # sigmas=np.divide(1,sigmacoeff)
    
    nvec=np.arange(N+1)
    coeff1=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    coeff2=np.multiply(qvec,qvec)
    coeff=np.multiply(coeff1,coeff2)
    const=dt/((M**3)*(a**4))
    sigmacoeff1=np.longdouble(2*const*coeff)*1000
    sigmacoeff=np.array(np.zeros(N+1),dtype=np.longdouble)
    
    sigmacoeff=1+sigmacoeff1
    print(sigmacoeff)
    #sigmacoeff=(1+2*p.dt*p.K4*coeff/p.a**4)
    #sigmas=np.array((np.divide(1,sigmacoeff)),dtype=np.longfloat)
    sigmas=np.array(((1/sigmacoeff)),dtype=np.longfloat)
    print(sigmas)
    
    for m in range(M+1):
        sigma[m,:]=sigmas
    
    
    
    return sigma

def sigmaPhiSV(M,N,qvec,a,dt):
    sigma=np.zeros((M+1,N+1))
    # nvec=np.arange(N+1)
    # coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
    # sigmacoeff=(1+2*dt*K4*coeff/a**4)
    # sigmas=np.divide(1,sigmacoeff)
    
    nvec=np.arange(N+1)
    coeff1=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))
    coeff2=np.multiply(qvec,qvec)
    coeff=np.multiply(coeff1,coeff2)
    const=dt/((M**3)*(a**4))
    
    sigmacoeff1=2*const*coeff
    sigmacoeff=(1+sigmacoeff1)
    
    #sigmacoeff=(1+2*p.dt*p.K4*coeff/p.a**4)
    sigmas=np.divide(1,sigmacoeff)
    
    for m in range(M+1):
        sigma[m,:]=sigmas
    
    print(sigmas)
    
    return sigma

def q(M,N):
    nc=int(M**(0.75)/1.5)
    
    qvec=np.zeros(N+1)
    for j in range(nc+1,N+1):
        qvec[j]=np.exp(-(j-M)**2/(2*(j-nc)**2))
    return qvec
