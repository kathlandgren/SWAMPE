# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:03:58 2020

@author: ek672
"""

import numpy as np
import scipy.special as sp

#local imports 
import params as p
#import explicit_time_diff as exp
import fft_legendre_trans as rfl
#import exp_t_diff_new as tdiff
import schwartztrauber as S

expflag=p.expflag
if expflag==1:
    import exp_t_diff_new as tdiff
else:
    import semi_imp_t_diff as tdiff


##TODO: fix inputs into invrsUV: fmn, tstepcoeffmn
def tstepping(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,marray,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi):

    newPhimn,newPhitstep2=tdiff.phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi)
    newdeltamn,newdeltatstep2=tdiff.delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi)
    newetamn,newetatstep2=tdiff.eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi)
    
    Unew,Vnew=rfl.invrsUV(newdeltamn,newetamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray) # TODO: reorder inputs later alphabetically
    
    return newetamn,newetatstep2,newdeltamn,newdeltatstep2,newPhimn,newPhitstep2,Unew,Vnew

def tstepping_latlon(delta0,delta1,eta0,eta1,Phi0,Phi1, w, mus):
    
    # 1 means "now", 0 means the previous time step
    
    U1,V1 = S.A20_A21(delta1,)
    
    
    return delta2, eta2, Ph2, U1, V1


def tstepcoeffmn(M,N,a):
    """
    Generates the coefficient that multiplies components in invrsUV
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
    tstepcoeff2=np.zeros((J,M+1))
    for m in range(M+1):
        tstepcoeff2[:,m]=(2*dt/a**2)
    return tstepcoeff2
    

def narray(M,N):
    narray=np.zeros((M+1,N+1))
    for m in range(M+1):
        for n in range(N+1):
            narray[m,n]=n*(n+1)
    return narray

def tstepcoeff(J,M,dt,mus,a):
    tstepcoeff=np.zeros((J,M+1))
    for m in range(M+1):
       tstepcoeff[:,m]=(2*dt)/(a*(1-mus**2))
    return tstepcoeff

def mJarray(J,M):
    mJarray=np.zeros((J,M+1))
    mtemp=np.arange(M+1)
    for j in range(J):
        mJarray[j,:]=mtemp
    return mJarray

def marray(M,N):
    marray=np.zeros((M+1,N+1)) #TODO make this an input, compute once
    mtemp=np.arange(M+1)
    for n in range(N+1):
        marray[:,n]=mtemp
    return marray