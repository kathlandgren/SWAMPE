# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 13:29:58 2021

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
def phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test):
    #use the "1" for state variables (the latest one)
    #tstepcoeff1=tstepcoeff1/2
    #tstepcoeff2=tstepcoeff2/2

    
    Phicomp1=rfl.fwd_leg(Phim0, J, M, N, Pmn, w)
    
    Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    Phicomp3prep=np.multiply(tstepcoeff1,Dm)
    Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
    Phicomp4=dt*Phibar*rfl.fwd_leg(deltam0, J, M, N, Pmn, w)
    
    
    deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
    

    deltacomp3prep=np.multiply(tstepcoeff1,Am)

    deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
    
    deltacomp5prep=np.multiply(tstepcoeff2,Em)
    deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    
    deltacomp5=np.multiply(narray,deltacomp5)
    
   
    Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))

    if forcflag==1:
        # Phiforcing=rfl.fwd_leg(2*dt*PhiFM, J, M, N, Pmn, w)
        # Phimntstep=Phimntstep+Phiforcing
        if test==10:
            deltaf3prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Fm)
            deltaf3=rfl.fwd_leg(deltaf3prep, J, M, N, Pmn, w)
            
            deltaf4prep=np.multiply(tstepcoeff1,Gm)
            deltaf4=rfl.fwd_leg(deltaf4prep, J, M, N, Hmn, w)
            
            deltaforcing=deltaf3+deltaf4
    
            Phiforcing=rfl.fwd_leg((dt)*PhiFM, J, M, N, Pmn, w)
    
            Phimntstep=Phimntstep+Phiforcing-(dt)*Phibar*deltaforcing   
            
         elif test==11:
            Phimntstep=Phimntstep+Phiforcing 
    
    if diffflag==1:
        Phimntstep=filters.diffusion(Phimntstep, sigmaPhi) 
        
        
    
    newPhimtstep=rfl.invrs_leg(Phimntstep, I,J, M, N, Pmn)
    newPhitstep=rfl.invrs_fft(newPhimtstep, I)

    
    return Phimntstep,newPhitstep

def delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test):
   # tstepcoeff1=tstepcoeff1/2
   # tstepcoeff2=tstepcoeff2/2
    
    
    deltacomp1=rfl.fwd_leg(deltam0, J, M, N, Pmn, w)
    
    deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
    

    deltacomp3prep=np.multiply(tstepcoeff1,Am)

    deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
    

    deltacomp4prep=np.multiply(tstepcoeff2,Phim0)
    deltacomp4=rfl.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
    
    deltacomp4=np.multiply(narray,deltacomp4)
    
    deltacomp5prep=np.multiply(tstepcoeff2,Em)
    deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    
    deltacomp5=np.multiply(narray,deltacomp5)
    
    
    #K1=deltacomp2+deltacomp3+deltacomp4+deltacomp5
    
    Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)


    Phicomp3prep=np.multiply(tstepcoeff1,Dm)
    Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)

    
    deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2

    if forcflag==1:
        if test==10:
        
            deltaf3prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Fm)
            deltaf3=rfl.fwd_leg(deltaf3prep, J, M, N, Pmn, w)
            
            deltaf4prep=np.multiply(tstepcoeff1,Gm)
            deltaf4=rfl.fwd_leg(deltaf4prep, J, M, N, Hmn, w)
            
            # deltaforcing=-deltaf1+deltaf2+deltaf3-deltaf4
            deltaforcing=deltaf3+deltaf4
    
            Phiforcing=np.multiply(narray,rfl.fwd_leg((dt)*PhiFM, J, M, N, Pmn, w))/a**2
    
            deltamntstep=deltamntstep+deltaforcing+Phiforcing
            

    
    if diffflag==1:
        deltamntstep=filters.diffusion(deltamntstep, sigma)

    newdeltamtstep=rfl.invrs_leg(deltamntstep, I,J, M, N, Pmn)
    newdeltatstep=rfl.invrs_fft(newdeltamtstep, I)
    return deltamntstep,newdeltatstep


def eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test):
    #tstepcoeff1=tstepcoeff1/2
    
    etacomp1=rfl.fwd_leg(etam0, J, M, N, Pmn, w)

    etacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Am)
    etacomp2=rfl.fwd_leg(etacomp2prep, J, M, N, Pmn, w)
    
    etacomp3prep=np.multiply(tstepcoeff1,Bm)
    etacomp3=rfl.fwd_leg(etacomp3prep, J, M, N, Hmn, w)
    
    
    
    etamntstep=etacomp1-etacomp2+etacomp3
    #etamntstep=etacomp1#-etacomp2+etacomp3
    
    if forcflag==1:
        
        if test==10:
            #eta forcing is explicit as in (A.53) in Hack and Jakob (1992)
            
            etaf3prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Gm)
            etaf3=rfl.fwd_leg(etaf3prep, J, M, N, Pmn, w)
            
            etaf4prep=np.multiply(tstepcoeff1,Fm)
            etaf4=rfl.fwd_leg(etaf4prep, J, M, N, Hmn, w)
            
    
            etaforcing=etaf3-etaf4
        
            etamntstep=etamntstep+etaforcing

    
    if diffflag==1:
        etamntstep=filters.diffusion(etamntstep, sigma)
    
    newetamtstep=rfl.invrs_leg(etamntstep, I,J, M, N, Pmn)
    newetatstep=rfl.invrs_fft(newetamtstep, I)
    return etamntstep,newetatstep
