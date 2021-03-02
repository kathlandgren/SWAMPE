# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:59:45 2020

@author: ek672
"""


import numpy as np
import matplotlib.pyplot as plt
#local imports
import params as p
import fft_legendre_trans as rfl
import filters
import testing_plots


## PHI tstep

def phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi):
    Phicoeff0=1/(1+(Phibar/4)*narray*(2*dt)**2/a**2)
    
    
    Phicoeff1=(1-(Phibar/4)*narray*(2*dt)**2/(a**2))
    Phicomp1prep=rfl.fwd_leg(Phim0, J, M, N, Pmn, w)
    
    Phicomp1=np.multiply(Phicoeff1,Phicomp1prep)
    
    Phicomp2prep=rfl.fwd_leg(deltam0,J,M,N,Pmn,w)
    Phicomp2=2*dt*Phibar*Phicomp2prep
    
    Phicomp3prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Pmn, w)
    
    Phicomp4prep=np.multiply(tstepcoeff1,Dm)
    Phicomp4=rfl.fwd_leg(Phicomp4prep, J, M, N, Hmn, w)
    
    Phicoeff5=narray*Phibar*(2*dt**2)/(2*a**2)
    Phicomp5prep=rfl.fwd_leg(Em,J,M,N,Pmn,w)
    Phicomp5=np.multiply(Phicoeff5,Phicomp5prep)
    
    Phicoeff6=Phibar*(2*dt)/2
    Phicomp6prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Bm))
    Phicomp6prep2=rfl.fwd_leg(Phicomp6prep, J, M, N, Pmn, w)
    Phicomp6=Phicoeff6*Phicomp6prep2
    
    Phicoeff7=Phicoeff6
    Phicomp7prep=np.multiply(tstepcoeff1,Am)
    Phicomp7prep2=rfl.fwd_leg(Phicomp7prep, J, M, N, Hmn, w)
    Phicomp7=Phicoeff7*Phicomp7prep2
    
    
    
    Phimntstep=Phicoeff0*(Phicomp1-Phicomp2-Phicomp3+Phicomp4-Phicomp5-Phicomp6-Phicomp7)

    if forcflag==1:
        Phiforcing=rfl.fwd_leg(2*dt*PhiFM, J, M, N, Pmn, w)
        Phimntstep=Phimntstep+Phiforcing
        
    
    if diffflag==1:
        Phimntstep=filters.diffusion(Phimntstep, sigmaPhi)    
    
    test,newPhimtstep=rfl.invrs_leg(Phimntstep, I,J, M, N, Pmn)
    newPhitstep=rfl.invrs_fft(newPhimtstep, I)

    
    return Phimntstep,newPhitstep

def delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi):
    deltacoeff1=(1-(Phibar/4)*narray*(2*dt)**2/(a**2))
    deltacomp1prep=rfl.fwd_leg(deltam0, J, M, N, Pmn, w)
    
    deltacomp1=np.multiply(deltacoeff1,deltacomp1prep)
    
    deltacoeff2=narray*2*dt/a**2
    deltacomp2prep=rfl.fwd_leg(Phim0+Em,J,M,M,Pmn,w)
    
    deltacomp2=np.multiply(deltacoeff2,deltacomp2prep)
    
    
    deltacomp3prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Pmn, w)
    

    deltacomp4prep=np.multiply(tstepcoeff1,Am)
    deltacomp4=rfl.fwd_leg(deltacomp4prep, J, M, N, Hmn, w)
    
    
    deltacomp5coeff=narray*dt/a**2
    deltacomp5prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Cm)
    deltacomp5prep2=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    deltacomp5=np.multiply(deltacomp5coeff,deltacomp5prep2)
    
    deltacomp6prep=np.multiply(tstepcoeff1,Dm)
    deltacomp6prep2=rfl.fwd_leg(deltacomp6prep, J, M, N, Hmn, w)
    deltacomp6=np.multiply(deltacomp5coeff,deltacomp6prep2)
    
    deltacoeff0=1/(1+(Phibar/4)*narray*(2*dt)**2/a**2)
    

    deltamntstep=deltacoeff0*(deltacomp1+deltacomp2+deltacomp3+deltacomp4-deltacomp5+deltacomp6)
    #deltamntstep=deltacoeff0*(deltacomp1)
    if forcflag==1:
        # if taudrag==-1:
        #     deltaf1=np.zeros((M+1,N+1))
        #     deltaf2=np.zeros((M+1,N+1))
        # else:
        #     deltaf1prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Um)/taudrag
        #     deltaf1=rfl.fwd_leg(deltaf1prep, J, M, N, Pmn, w)
        
        #     deltaf2prep=np.multiply(tstepcoeff1,Vm)/taudrag
        #     deltaf2=rfl.fwd_leg(deltaf2prep, J, M, N, Hmn, w)
       
        
        deltaf3prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Fm)
        deltaf3=rfl.fwd_leg(deltaf3prep, J, M, N, Pmn, w)
        
        deltaf4prep=np.multiply(tstepcoeff1,Gm)
        deltaf4=rfl.fwd_leg(deltaf4prep, J, M, N, Hmn, w)
        
        # deltaforcing=-deltaf1+deltaf2+deltaf3-deltaf4
        deltaforcing=deltaf3+deltaf4
        deltamntstep=deltamntstep+deltaforcing
        
    
    if diffflag==1:
        deltamntstep=filters.diffusion(deltamntstep, sigma)

    test,newdeltamtstep=rfl.invrs_leg(deltamntstep, I,J, M, N, Pmn)
    newdeltatstep=rfl.invrs_fft(newdeltamtstep, I)
    return deltamntstep,newdeltatstep


def eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi):

    ## ETA tstep seems to work!
    etacomp1=rfl.fwd_leg(etam0, J, M, N, Pmn, w)

    etacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Am)
    etacomp2=rfl.fwd_leg(etacomp2prep, J, M, N, Pmn, w)
    
    etacomp3prep=np.multiply(tstepcoeff1,Bm)
    etacomp3=rfl.fwd_leg(etacomp3prep, J, M, N, Hmn, w)
    
    
    
    etamntstep=etacomp1-etacomp2+etacomp3
    #etamntstep=etacomp1#-etacomp2+etacomp3
    
    if forcflag==1:
        # if taudrag==-1:
        #     etaf1=np.zeros((M+1,N+1))
        #     etaf2=np.zeros((M+1,N+1))
        # else:
        #     etaf1prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Vm)/taudrag
        #     etaf1=rfl.fwd_leg(etaf1prep, J, M, N, Pmn, w)
            
        #     etaf2prep=np.multiply(tstepcoeff1,Um)/taudrag
        #     etaf2=rfl.fwd_leg(etaf2prep, J, M, N, Hmn, w)
        
        etaf3prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Gm)
        etaf3=rfl.fwd_leg(etaf3prep, J, M, N, Pmn, w)
        
        etaf4prep=np.multiply(tstepcoeff1,Fm)
        etaf4=rfl.fwd_leg(etaf4prep, J, M, N, Hmn, w)
        
        #etaforcing=-etaf1+etaf2+etaf3+etaf4
        etaforcing=etaf3-etaf4
    
        etamntstep=etamntstep+etaforcing

    
    if diffflag==1:
        etamntstep=filters.diffusion(etamntstep, sigma)
    
    test,newetamtstep=rfl.invrs_leg(etamntstep, I,J, M, N, Pmn)
    newetatstep=rfl.invrs_fft(newetamtstep, I)
    return etamntstep,newetatstep
