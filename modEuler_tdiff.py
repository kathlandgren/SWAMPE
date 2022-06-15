# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:48:34 2021

@author: ek672
"""


import numpy as np
import scipy.special as sp
#local imports
import spectral_transform as rfl
import filters


## PHI tstep
def phi_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
    #use the "1" for state variables (the latest one)
    tstepcoeff1=tstepcoeff1/2
    tstepcoeff2=tstepcoeff2/2

    
    Phicomp1=rfl.fwd_leg(Phim1, J, M, N, Pmn, w)
    
    Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    Phicomp3prep=np.multiply(tstepcoeff1,Dm)
    Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
    Phicomp4=dt*Phibar*rfl.fwd_leg(deltam1, J, M, N, Pmn, w)
    
    
    deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
    

    deltacomp3prep=np.multiply(tstepcoeff1,Am)

    deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
    
    deltacomp5prep=np.multiply(tstepcoeff2,Em)
    deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    
    deltacomp5=np.multiply(narray,deltacomp5)
    # print('min Phicomp2 '+str(np.min(Phicomp2)))
    # print('max Phicomp3 '+str(np.max(Phicomp3)))
    # print('min Phicomp4 '+str(np.min(Phicomp4)))
    
   
    Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))

    if forcflag==1:
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
    
        
        Phicomp1=rfl.fwd_leg(Phim1, J, M, N, Pmn, w)
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
        
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
        
        Phicomp4=dt*Phibar*rfl.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm+Fm)
        deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am-Gm)
    
        deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
        # print('min Phicomp2 '+str(np.min(Phicomp2)))
        # print('max Phicomp3 '+str(np.max(Phicomp3)))
        # print('min Phicomp4 '+str(np.min(Phicomp4)))
        
       
        Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))


        Phiforcing=rfl.fwd_leg((dt)*PhiFM, J, M, N, Pmn, w)

        Phimntstep=Phimntstep+Phiforcing
    else:
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
    
        
        Phicomp1=rfl.fwd_leg(Phim1, J, M, N, Pmn, w)
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
        
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
        
        Phicomp4=dt*Phibar*rfl.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
        deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am)
    
        deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
        # print('min Phicomp2 '+str(np.min(Phicomp2)))
        # print('max Phicomp3 '+str(np.max(Phicomp3)))
        # print('min Phicomp4 '+str(np.min(Phicomp4)))
        
       
        Phimntstep=Phicomp1-Phicomp2+Phicomp3-Phicomp4-Phibar*(0.5)*(deltacomp2+deltacomp3+deltacomp5+(1/a**2)*np.multiply(narray,Phicomp1))


    
    if diffflag==1:
        Phimntstep=filters.diffusion(Phimntstep, sigmaPhi) 
        
        
    
    newPhimtstep=rfl.invrs_leg(Phimntstep, I,J, M, N, Pmn)
    newPhitstep=rfl.invrs_fft(newPhimtstep, I)

    
    return Phimntstep,newPhitstep

def delta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
    tstepcoeff1=tstepcoeff1/2
    tstepcoeff2=tstepcoeff2/2
    
    
    deltacomp1=rfl.fwd_leg(deltam1, J, M, N, Pmn, w)
    
    deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm)
    deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
    

    deltacomp3prep=np.multiply(tstepcoeff1,Am)

    deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
    
    if test==9:
        Phis=np.zeros((J,I))
        [mus,w]=sp.roots_legendre(J)
        for i in range(I):
            for j in range(J):
                Phis[j,:]=1500*9.8*(1-(mus[j])**2)
       

        if t*dt<1*24*3600:
            factor=t*dt/(1*24*3600)
        else:
            factor=1    
        Phis=factor*Phis
            
        Phism=rfl.fwd_fft_trunc(Phis, I, M) 
        deltacomp4prep=np.multiply(tstepcoeff2,Phim1+Phism)
    else:
        deltacomp4prep=np.multiply(tstepcoeff2,Phim1)
    deltacomp4=rfl.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
    
    deltacomp4=np.multiply(narray,deltacomp4)
    
    deltacomp5prep=np.multiply(tstepcoeff2,Em)
    deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
    
    deltacomp5=np.multiply(narray,deltacomp5)
    

    
    Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
    Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)


    Phicomp3prep=np.multiply(tstepcoeff1,Dm)
    Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)

    # print('max deltacomp2 '+str(np.max(deltacomp2)))
    # print('max deltacomp3 '+str(np.max(deltacomp3)))
    # print('max deltacomp4 '+str(np.max(deltacomp4)))
    # print('max deltacomp5 '+str(np.max(deltacomp5)))
    
    # print('max remaining deltacomp '+str(np.max(np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2)))
    
    deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2

    # print('max deltamn tstep '+str(np.max(deltamntstep)))
    if forcflag==1:
        
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
        
        
        deltacomp1=rfl.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm+Fm)
        deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am-Gm)
    
        deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        if test==9:
            Phis=np.zeros((J,I))
            [mus,w]=sp.roots_legendre(J)
            for i in range(I):
                for j in range(J):
                    Phis[j,:]=1500*9.8*(1-(mus[j])**2)
           
    
            if t*dt<1*24*3600:
                factor=t*dt/(1*24*3600)
            else:
                factor=1    
            Phis=factor*Phis
                
            Phism=rfl.fwd_fft_trunc(Phis, I, M) 
            deltacomp4prep=np.multiply(tstepcoeff2,Phim1+Phism)
        else:
            deltacomp4prep=np.multiply(tstepcoeff2,Phim1)
        deltacomp4=rfl.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
        
        deltacomp4=np.multiply(narray,deltacomp4)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
        
    
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
        # print('max deltacomp2 '+str(np.max(deltacomp2)))
        # print('max deltacomp3 '+str(np.max(deltacomp3)))
        # print('max deltacomp4 '+str(np.max(deltacomp4)))
        # print('max deltacomp5 '+str(np.max(deltacomp5)))
        
        # print('max remaining deltacomp '+str(np.max(np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2)))
        
        deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2

    # print('max deltamn tstep '+str(np.max(deltamntstep)))
        Phiforcing=np.multiply(narray,rfl.fwd_leg((dt/2)*PhiFM, J, M, N, Pmn, w))/a**2

        deltamntstep=deltamntstep+Phiforcing
        
    else:
            
        
        tstepcoeff1=tstepcoeff1/2
        tstepcoeff2=tstepcoeff2/2
        
        
        deltacomp1=rfl.fwd_leg(deltam1, J, M, N, Pmn, w)
        
        deltacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Bm+Fm)
        deltacomp2=rfl.fwd_leg(deltacomp2prep, J, M, N, Pmn, w)
        
    
        deltacomp3prep=np.multiply(tstepcoeff1,Am-Gm)
    
        deltacomp3=rfl.fwd_leg(deltacomp3prep, J, M, N, Hmn, w)
        
        if test==9:
            Phis=np.zeros((J,I))
            [mus,w]=sp.roots_legendre(J)
            for i in range(I):
                for j in range(J):
                    Phis[j,:]=1500*9.8*(1-(mus[j])**2)
           
    
            if t*dt<1*24*3600:
                factor=t*dt/(1*24*3600)
            else:
                factor=1    
            Phis=factor*Phis
                
            Phism=rfl.fwd_fft_trunc(Phis, I, M) 
            deltacomp4prep=np.multiply(tstepcoeff2,Phim1+Phism)
        else:
            deltacomp4prep=np.multiply(tstepcoeff2,Phim1)
        deltacomp4=rfl.fwd_leg(deltacomp4prep, J, M, N, Pmn, w)
        
        deltacomp4=np.multiply(narray,deltacomp4)
        
        deltacomp5prep=np.multiply(tstepcoeff2,Em)
        deltacomp5=rfl.fwd_leg(deltacomp5prep, J, M, N, Pmn, w)
        
        deltacomp5=np.multiply(narray,deltacomp5)
        
    
        
        Phicomp2prep=np.multiply(tstepcoeff1,np.multiply((1j)*mJarray,Cm))
        Phicomp2=rfl.fwd_leg(Phicomp2prep, J, M, N, Pmn, w)
    
    
        Phicomp3prep=np.multiply(tstepcoeff1,Dm)
        Phicomp3=rfl.fwd_leg(Phicomp3prep, J, M, N, Hmn, w)
    
        # print('max deltacomp2 '+str(np.max(deltacomp2)))
        # print('max deltacomp3 '+str(np.max(deltacomp3)))
        # print('max deltacomp4 '+str(np.max(deltacomp4)))
        # print('max deltacomp5 '+str(np.max(deltacomp5)))
        
        # print('max remaining deltacomp '+str(np.max(np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2)))
        
        deltamntstep=deltacomp1+deltacomp2+deltacomp3+deltacomp4+deltacomp5+np.multiply(narray,(Phicomp2+Phicomp3)/2)/a**2-Phibar*np.multiply(narray,deltacomp1)/a**2
   
 
        
        

    
    if diffflag==1:
        deltamntstep=filters.diffusion(deltamntstep, sigma)

    newdeltamtstep=rfl.invrs_leg(deltamntstep, I,J, M, N, Pmn)
    newdeltatstep=rfl.invrs_fft(newdeltamtstep, I)
    return deltamntstep,newdeltatstep


def eta_timestep(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,Pmn,Hmn,w,tstepcoeff1,tstepcoeff2,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi,test,t):
 
    
    if forcflag==1:
        etacomp1=rfl.fwd_leg(etam1, J, M, N, Pmn, w)
    
        etacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Am-Gm)
        etacomp2=rfl.fwd_leg(etacomp2prep, J, M, N, Pmn, w)
        
        etacomp3prep=np.multiply(tstepcoeff1,Bm+Fm)
        etacomp3=rfl.fwd_leg(etacomp3prep, J, M, N, Hmn, w)
        etamntstep=etacomp1-etacomp2+etacomp3
        
    else:
        tstepcoeff1=tstepcoeff1/2

        etacomp1=rfl.fwd_leg(etam1, J, M, N, Pmn, w)
    
        etacomp2prep=np.multiply(np.multiply(tstepcoeff1,(1j)*mJarray),Am)
        etacomp2=rfl.fwd_leg(etacomp2prep, J, M, N, Pmn, w)
        
        etacomp3prep=np.multiply(tstepcoeff1,Bm)
        etacomp3=rfl.fwd_leg(etacomp3prep, J, M, N, Hmn, w)
    
    
    
    etamntstep=etacomp1-etacomp2+etacomp3
    #etamntstep=etacomp1#-etacomp2+etacomp3
    
    
   
    
    if diffflag==1:
        etamntstep=filters.diffusion(etamntstep, sigma)
    
    newetamtstep=rfl.invrs_leg(etamntstep, I,J, M, N, Pmn)
    newetatstep=rfl.invrs_fft(newetamtstep, I)
    return etamntstep,newetatstep
