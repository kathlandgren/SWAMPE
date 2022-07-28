# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:48:34 2021

@author: ek672
"""


import numpy as np
import scipy.special as sp
#local imports
from . import spectral_transform as st
from . import filters


## PHI tstep
def phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
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

def delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
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


def eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFm,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
 
    
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
