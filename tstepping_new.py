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
import pyshtools as pysh
import forcing as f

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

def tstepping_latlon(test,U0,V0,delta0,delta1,zeta0,zeta1,f_latlon,Phi0,Phi1, w, mus,J,M,nMAT1,nMAT2,nMAT3,mnMAT1,mnMAT2,mnMAT3,mnMAT4,mnMAT5,musMAT,a,dt,Phibar,normnum,forcflag,PhiF,F,G):
    
    # 1 means "now", 0 means the previous time step
    
    #7.1: fwrd transform zeta and delta    
    #7.2: get U, V from the above
    if test<2:
        U1=U0
        V1=V0
    else: #elif test==10:
        U1,V1=S.A20_A21(delta1,zeta1,M,nMAT3,mnMAT1,mnMAT2,mnMAT3,w,mus,J,normnum)
    #7.3: make zeta and delta and forward transfrom, get RHS
    X=np.multiply(zeta1+f_latlon,V1) #are we setting up the vorticity twice here? 
    Y=np.multiply(-(zeta1+f_latlon),U1)
    
    brmn, bimn=S.A22_A23(X,Y,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)
    crmn, cimn=S.A24_A25(X,Y,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)
    
    zetaRHS=S.A14(crmn,cimn,nMAT1,mus,M,J,normnum)
    deltaRHS1=S.A15(brmn,bimn,nMAT1,mus,M,J,normnum)
    
    #7.4 get the geopotential RHS
    X=np.multiply(Phi1,U1)
    Y=np.multiply(Phi1,V1)
    
    brmn, bimn=S.A22_A23(X,Y,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)
    
    PhiRHS=S.A15(brmn,bimn,nMAT1,mus,M,J,normnum)
    
    
    #7.5
    
    deltaRHS2=S.step7p5(Phi1,U1,V1,w,mus,J,M,musMAT,nMAT2,normnum)
    
    #make outermost coefficient 
    acosMAT=(1/a)*musMAT
    #print(np.shape(musMAT))
    #print(np.shape(zetaRHS))
    #timestepping
    zeta2=zeta0 +(2*dt)*(-np.multiply(acosMAT,zetaRHS))
    delta2=delta0+(2*dt)*(np.multiply(acosMAT,deltaRHS1)-deltaRHS2/a**2)
    Phi2=Phi0+(2*dt)*(-np.multiply(acosMAT,PhiRHS))
    #Phi2=Phi0+(2*dt)*(-np.multiply(acosMAT,PhiRHS+Phibar*delta1))
    if forcflag==1:
        Phi2=Phi2+(2*dt)*PhiF#*(np.multiply(acosMAT,PhiF))
        
        Fbrmn, Fbimn=S.A22_A23(F,G,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)
        Fcrmn, Fcimn=S.A24_A25(F,G,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)
        zetaRHSF=S.A14(Fcrmn,Fcimn,nMAT1,mus,M,J,normnum)  
        deltaRHSF=S.A15(Fbrmn,Fbimn,nMAT1,mus,M,J,normnum)
        
        delta2=delta2+(2*dt)*(np.multiply(acosMAT,deltaRHSF))
        zeta2=zeta2+(2*dt)*(np.multiply(acosMAT,zetaRHSF))
        
    print(np.max(PhiRHS))
    print(np.max(zetaRHS))
    print(np.max(deltaRHS1 - deltaRHS2/a))
    
    # # For 1/cos^2(lat) in the integrals A22-25
    # # make outermost coefficient 
    # leapfrog=2*dt/a

    # #timestepping
    # zeta2=zeta0 - (leapfrog)*(zetaRHS)
    # delta2=delta0 + (leapfrog)*(deltaRHS1 - deltaRHS2/a)
    # Phi2=Phi0 - (leapfrog)*(PhiRHS)
    # #Phi2=Phi0-(leapfrog)*(PhiRHS + a*Phibar*delta1)
    
    #filter zeta1, delta1, Phi1
    
    return delta2, zeta2, Phi2, U1, V1


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