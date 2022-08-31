"""
This module contains the main SWAMPE function which runs the 2D shallow-water general circulation model.
"""
# Import python packages
import numpy as np
# Import program packages
from . import initial_conditions 
from . import spectral_transform 
from . import time_stepping 
from . import plotting 
from . import forcing 
from . import filters 
from . import continuation

def run_model(M,dt,tmax,Phibar, omega, a, test=None, g=9.8, forcflag=True, taurad=86400, taudrag=86400, DPhieq=4*(10**6), a1=0.05, plotflag=True, plotfreq=5, minlevel=None, maxlevel=None, diffflag=True,modalflag=True,alpha=0.01,contflag=False,saveflag=True,expflag=False,savefreq=150, K6=1.24*10**33,custompath=None,contTime=None,timeunits='hours',verbose=True):    
    """_summary_

    :param M: spectral resolution
    :type M: int
    :param dt: time step length in seconds
    :type dt: float
    :param tmax: number of timesteps to run in the simulation
    :type tmax: int
    :param Phibar: reference geopotential (a good rule of thumb is Phibar=gH, where H is scale height in meters)
    :type Phibar: float
    :param omega: planetary rotation rate in radians/s
    :type omega: float
    :param a: planetary radius in meters
    :type a: float
    :param test: number of test from Jakob & Hack (1993), tests 1 and 2 are supported, defaults to None
    :type test: int, optional
    :param g: surface gravity, defaults to 9.8 m/s^2
    :type g: float, optional
    :param forcflag: option to implement radiative forcing from the host star, defaults to True
    :type forcflag: bool, optional
    :param taurad: radiative timescale for Newtonian relaxation in the forcing scheme, defaults to 86400 s
    :type taurad: float, optional
    :param taudrag: drag timescale for wind forcing, defaults to 86400 s
    :type taudrag: float, optional
    :param DPhieq: day/night amplitude of prescribed local radiative equilibrium geopotential, defaults to 4*(10**6) m^2/s^2.
    :type DPhieq: float, optional
    :param a1: angle for tests 1 and 2 from Jakob and Hack (1993), defaults to 0.05
    :type a1: float, optional
    :param plotflag: option to display progress plots over the course of the simulation run, defaults to True
    :type plotflag: bool, optional
    :param plotfreq: frequency of plot output during the simulation run, defaults to once every 5 timesteps
    :type plotfreq: int, optional
    :param minlevel: minimum level of colorbar for geopotential plotting, defaults to minimum geopotential
    :type minlevel: float, optional
    :param maxlevel: maximum level of colorbar for geopotential plotting, defaults to minimum geopotential
    :type maxlevel: float, optional
    :param diffflag: option to turn on the diffusion/hyperviscosity filter, defaults to True (strongly recommended)
    :type diffflag: bool, optional
    :param modalflag: option to turn on the modal splitting filter from Hack and Jakob (1992), defaults to True
    :type modalflag: bool, optional
    :param alpha: parameter for the modal splitting filter from Hack and Jakob (1992), defaults to 0.01
    :type alpha: float, optional
    :param contflag: option to continue the simulation from saved data, defaults to False
    :type contflag: bool, optional
    :param saveflag: option to save data as pickle files, defaults to True
    :type saveflag: bool, optional
    :param expflag: option to use explicit time-stepping scheme instead of modified Euler, defaults to False (strongly recommended to use the modified Euler scheme)
    :type expflag: bool, optional
    :param savefreq: frequency of saving the data, defaults to once every 150 timesteps
    :type savefreq: int, optional
    :param K6: sixth order hyperviscosity filter parameter, defaults to 1.24*10**33
    :type K6: float, optional
    :param custompath: option to specify the path for saving the data. By default SWAMPE will make a local directory data/ and store files there.
    :type custompath: str, optional
    :param contTime: (if continuing from saved data) timestamp of the data, defaults to None
    :type contTime: int, optional
    :param timeunits: time units, defaults to 'hours', also supports 'minutes' and 'seconds'
    :type timeunits: str, optional
    :param verbose: option to print progress statements, defaults to True
    :type verbose: bool, optional
    """

    #get other dimensional parameters using the spectral dimension
    N,I,J,otherdt,lambdas,mus,w=initial_conditions.spectral_params(M)
    
    

    # Associated Legendre Polynomials and their derivatives
    Pmn, Hmn = spectral_transform.PmnHmn(J, M, N, mus)
    
    sigma=filters.sigma6(M,N,K6,a, dt)
    sigmaPhi=filters.sigma6Phi(M, N, K6, a, dt)
        
    ## Initialize data arrays 
   
    etadata=np.zeros((3,J,I))
    deltadata=np.zeros((3,J,I))
    Phidata=np.zeros((3,J,I))
    
    etamdata=np.zeros((3,J,M+1),dtype=complex)
    deltamdata=np.zeros((3,J,M+1),dtype=complex)
    Phimdata=np.zeros((3,J,M+1),dtype=complex)
    
    etamndata=np.zeros((3,M+1,N+1),dtype=complex)
    deltamndata=np.zeros((3,M+1,N+1),dtype=complex)
    Phimndata=np.zeros((3,M+1,N+1),dtype=complex)
    
    Udata=np.zeros((3,J,I))
    Vdata=np.zeros((3,J,I))
    
    Umdata=np.zeros((3,J,M+1),dtype=complex)
    Vmdata=np.zeros((3,J,M+1),dtype=complex)
    
    Adata=np.zeros((3,J,I))
    Bdata=np.zeros((3,J,I))
    Cdata=np.zeros((3,J,I))
    Ddata=np.zeros((3,J,I))
    Edata=np.zeros((3,J,I))
    
    Amdata=np.zeros((3,J,M+1),dtype=complex)
    Bmdata=np.zeros((3,J,M+1),dtype=complex)
    Cmdata=np.zeros((3,J,M+1),dtype=complex)
    Dmdata=np.zeros((3,J,M+1),dtype=complex)
    Emdata=np.zeros((3,J,M+1),dtype=complex)
    
    Fdata=np.zeros((3,J,I))
    Gdata=np.zeros((3,J,I))
    
    Fmdata=np.zeros((3,J,M+1),dtype=complex)
    Gmdata=np.zeros((3,J,M+1),dtype=complex)
    
    Phiforcingdata=np.zeros((3,J,I))
    Phiforcingmdata=np.zeros((3,J,M+1),dtype=complex)
    
    #long arrays
    spinupdata=np.zeros((tmax,2))
    
    geopotdata=np.zeros((tmax,2))
    
    
    #coriolis
    fmn=initial_conditions.coriolismn(M, omega)
    
    tstepcoeffmn=time_stepping.tstepcoeffmn(M,N,a)
    tstepcoeff=time_stepping.tstepcoeff(J,M,dt,mus,a)
    tstepcoeff2=time_stepping.tstepcoeff2(J,M,dt,a)
    mJarray=time_stepping.mJarray(J,M)
    marray=time_stepping.marray(M, N)
    narray=time_stepping.narray(M,N)
        
    

    if contflag==False:
        SU0, sina, cosa, etaamp,Phiamp=initial_conditions.test1_init(a, omega, a1)
        
        if test==1:
            etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=initial_conditions.state_var_init(I,J,mus,lambdas,test,etaamp,a,sina,cosa,Phibar,Phiamp)
        elif test==2:
            etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=initial_conditions.state_var_init(I,J,mus,lambdas,test,etaamp,a,sina,cosa,Phibar,Phiamp)   
        else:
            etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=initial_conditions.state_var_init(I,J,mus,lambdas,test,etaamp)
        
        Uic,Vic=initial_conditions.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
    
    elif contflag==True:
        
        if custompath==None:

            etaic0=continuation.read_pickle('eta-'+str(contTime))
            etaic1 = etaic0

            deltaic0=continuation.read_pickle('delta-'+str(contTime))
            deltaic1 = deltaic0

            Phiic0=continuation.read_pickle('Phi-'+str(contTime))
            Phiic1 = Phiic0
        else:
            etaic0=continuation.read_pickle('eta-'+str(contTime),custompath=custompath)
            etaic1 = etaic0

            deltaic0=continuation.read_pickle('delta-'+str(contTime),custompath=custompath)
            deltaic1 = deltaic0

            Phiic0=continuation.read_pickle('Phi-'+str(contTime),custompath=custompath)
            Phiic1 = Phiic0
            

        
        etam0=spectral_transform.fwd_fft_trunc(etaic0, I, M)
        etamn0=spectral_transform.fwd_leg(etam0,J,M,N,Pmn,w)
        deltam0=spectral_transform.fwd_fft_trunc(deltaic0, I, M)
        deltamn0=spectral_transform.fwd_leg(deltam0,J,M,N,Pmn,w)
        
        Uiccomp,Viccomp=spectral_transform.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
        Uic=np.real(Uiccomp)
        Vic=np.real(Viccomp)
        
    Aic,Bic,Cic,Dic,Eic=initial_conditions.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)
    
    
    ## Store initial conditions in the data arrays files for easy access
    etadata[0,:,:]=etaic0
    etadata[1,:,:]=etaic1
    
    deltadata[0,:,:]=deltaic0
    deltadata[1,:,:]=deltaic1
    
    Phidata[0,:,:]=Phiic0
    Phidata[1,:,:]=Phiic1
    
    Udata[0,:,:]=Uic
    Udata[1,:,:]=Uic
    
    Vdata[0,:,:]=Vic
    Vdata[1,:,:]=Vic
    
    Adata[0,:,:]=Aic
    Adata[1,:,:]=Aic
    
    Bdata[0,:,:]=Bic
    Bdata[1,:,:]=Bic
    
    Cdata[0,:,:]=Cic
    Cdata[1,:,:]=Cic
    
    Ddata[0,:,:]=Dic
    Ddata[1,:,:]=Dic
    
    Edata[0,:,:]=Eic
    Edata[1,:,:]=Eic
    
    
    
    # Spin-up RMS wind calculations
    if contflag==False:
        spinupdata[0,0] = np.min(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))   
        spinupdata[0,1]=time_stepping.RMS_winds(a, I, J, lambdas, mus, Udata[0,:,:], Vdata[0,:,:])
        
        
        #geopotential calculations
        
        geopotdata[0,0]=np.min(Phidata[0,:,:])
        geopotdata[0,1]=np.max(Phidata[0,:,:])
    
    
    #### Forcing ####
    
    if test==None:
        Phieq=forcing.Phieqfun(Phibar, DPhieq, lambdas, mus, I, J, g)
        
        Q=forcing.Qfun(Phieq, Phiic0, Phibar,taurad)
        
        #geopotential forcing to be passed to time stepping
        PhiF=Q
        
        Phiforcingdata[0,:,:]=PhiF
        Phiforcingdata[1,:,:]=PhiF
        
        
        F,G=forcing.Rfun(Uic, Vic, Q, Phiic0,Phibar,taudrag)
        Fdata[0,:,:]=F
        Fdata[1,:,:]=Fdata[0,:,:]
        Gdata[0,:,:]=G
        Gdata[1,:,:]=Gdata[0,:,:]
 
    
    
    ####
    #Forward Fourier and Legendre transform
    ####
    
    ## FFT    
    
    Amdata[0,:,:]=spectral_transform.fwd_fft_trunc(Aic, I, M)
    Amdata[1,:,:]=Amdata[0,:,:]
    
    Bmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Bic, I, M)
    Bmdata[1,:,:]=Bmdata[0,:,:]
    
    Cmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Cic, I, M)
    Cmdata[1,:,:]=Cmdata[0,:,:]
    
    Dmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Dic, I, M)
    Dmdata[1,:,:]=Dmdata[0,:,:]
    
    Emdata[0,:,:]=spectral_transform.fwd_fft_trunc(Eic, I, M)
    Emdata[1,:,:]=Emdata[0,:,:]
    
    etamdata[0,:,:]=spectral_transform.fwd_fft_trunc(etaic0, I, M)
    etamdata[1,:,:]=etamdata[0,:,:]
    
    deltamdata[0,:,:]=spectral_transform.fwd_fft_trunc(deltaic0, I, M)
    deltamdata[1,:,:]=deltamdata[0,:,:]
    
    Phimdata[0,:,:]=spectral_transform.fwd_fft_trunc(Phiic0, I, M)
    Phimdata[1,:,:]=Phimdata[0,:,:]
    
    
    ## Forcing Fourier transform ##
    
    Phiforcingmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Phiforcingdata[0,:,:], I, M)
    Phiforcingmdata[1,:,:]=Phiforcingmdata[0,:,:]
    
    Fmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Fdata[0,:,:], I, M)
    Fmdata[1,:,:]=Fmdata[0,:,:]
    
    Gmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Gdata[0,:,:], I, M)
    Gmdata[1,:,:]=Gmdata[0,:,:]
    
    Umdata[0,:,:]=spectral_transform.fwd_fft_trunc(Udata[0,:,:], I, M)
    Umdata[1,:,:]=Umdata[0,:,:]
    
    Vmdata[0,:,:]=spectral_transform.fwd_fft_trunc(Vdata[0,:,:], I, M)
    Vmdata[1,:,:]=Vmdata[0,:,:]
    
    
    
    ## Forward Legendre
    
    etamndata[0,:,:]=spectral_transform.fwd_leg(etamdata[0,:,:],J,M,N,Pmn,w)
    etamndata[1,:,:]=etamndata[0,:,:]
    
    deltamndata[0,:,:]=spectral_transform.fwd_leg(deltamdata[0,:,:],J,M,N,Pmn,w)
    deltamndata[1,:,:]=deltamndata[0,:,:]
    
    Phimndata[0,:,:]=spectral_transform.fwd_leg(Phimdata[0,:,:],J,M,N,Pmn,w)
    Phimndata[1,:,:]=Phimndata[0,:,:]
    
    ####
    # Time stepping
    ####
    
    if contflag==False:
        starttime=2
    elif contflag==True:
        starttime=continuation.compute_t_from_timestamp(timeunits,contTime,dt)
        
    
    for t in range(starttime,tmax):
        
        etam0=etamdata[0,:,:]
        etam1=etamdata[1,:,:]
        
        deltam0=deltamdata[0,:,:]
        deltam1=deltamdata[1,:,:]
        
        Phim0=Phimdata[0,:,:]
        Phim1=Phimdata[1,:,:]
        
        Am=Amdata[1,:,:]
        Bm=Bmdata[1,:,:]
        Cm=Cmdata[1,:,:]    
        Dm=Dmdata[1,:,:]
        Em=Emdata[1,:,:]
        
        
        #forcing 
        Um=Umdata[1,:,:]
        Vm=Vmdata[1,:,:]
        
        Fm=Fmdata[1,:,:]
        Gm=Gmdata[1,:,:]
        
    
        PhiFm=Phiforcingmdata[1,:,:]    
        
        newetamn,neweta,newdeltamn,newdelta,newPhimn,newPhi,newU,newV=time_stepping.tstepping(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,marray,mJarray,narray,PhiFm,dt,a,Phibar,taurad,taudrag,forcflag,diffflag,expflag,sigma,sigmaPhi,test,t)
        
        
        #write new data
        etamndata[2,:,:]=newetamn
        deltamndata[2,:,:]=newdeltamn
        Phimndata[2,:,:]=newPhimn
            
        etadata[2,:,:]=np.real(neweta)
        deltadata[2,:,:]=np.real(newdelta)
        Phidata[2,:,:]=np.real(newPhi)
        
            
        if modalflag==True:
            if t>2:

                temp=np.zeros((J,I))
                temp[:,:]=Phidata[1,:,:]
                Phidata[1,:,:]=filters.modal_splitting(Phidata,alpha)
                etadata[1,:,:]=filters.modal_splitting(etadata,alpha)
                deltadata[1,:,:]=filters.modal_splitting(deltadata,alpha)


        
    
        if test==1:
            newU=Uic
            newV=Vic
    
    
        
        Udata[2,:,:]=np.real(newU)
        Vdata[2,:,:]=np.real(newV)
        
        spinupdata[t-1,0] = np.min(np.sqrt(Udata[1,:,:]**2 + Vdata[1,:,:]**2 ))
        spinupdata[t-1,1] = time_stepping.RMS_winds(a, I, J, lambdas, mus, Udata[1,:,:], Vdata[1,:,:])
        
        
        geopotdata[t-1,0]=np.min(Phidata[1,:,:])
        geopotdata[t-1,1]=np.max(Phidata[1,:,:])
        
        
        if saveflag==True:
   
            if dt*t%savefreq==0:                
                    
                timestamp=continuation.compute_timestamp(timeunits,dt,t)
                
                continuation.save_data(timestamp,np.real(neweta),np.real(newdelta), np.real(newPhi), np.real(newU), np.real(newV), spinupdata,geopotdata,custompath=custompath)

        
        if spinupdata[t-1,1]>8000:
            print('Time stepping stopped due to wind blow up. Max RMS winds = '+str(spinupdata[t-1,1]))
            break
            
            #t=tmax
    
        
        Umdata[2,:,:]=spectral_transform.fwd_fft_trunc(newU,I,M)
        Vmdata[2,:,:]=spectral_transform.fwd_fft_trunc(newV,I, M)
        
        Um=Umdata[2,:,:]
        Vm=Vmdata[2,:,:]
        
        
        etamdata[2,:,:]=spectral_transform.fwd_fft_trunc(neweta,I,M)
        deltamdata[2,:,:]=spectral_transform.fwd_fft_trunc(newdelta,I,M)
        Phimdata[2,:,:]=spectral_transform.fwd_fft_trunc(newPhi,I,M)
        
            
        ######## FORCING ############
        
        if test==None:
            
            Q=forcing.Qfun(Phieq, np.real(newPhi),Phibar, taurad)
            #geopotential forcing to be passed to time stepping
            PhiF=Q

            F,G=forcing.Rfun(np.real(newU), np.real(newV), Q, np.real(newPhi),Phibar,taudrag)
            
            Fdata[2,:,:]=F
            Gdata[2,:,:]=G
            
            Fmdata[2,:,:]=spectral_transform.fwd_fft_trunc(F, I, M)
            Gmdata[2,:,:]=spectral_transform.fwd_fft_trunc(G, I, M)
            
            Phiforcingdata[2,:,:]=PhiF
            Phiforcingmdata[2,:,:]=spectral_transform.fwd_fft_trunc(PhiF, I, M)  
        
            
            
        if verbose==True:
            if t%int(tmax/10)==0:
                print('t='+str(t)+', '+str(t*100/tmax)+'% complete')
        
        if plotflag==True:
            
            if t%plotfreq==0:

                
                timestamp=continuation.compute_timestamp(timeunits,dt,t)
                fig_zonal=plotting.mean_zonal_wind_plot(Udata[2,:,:], mus, timestamp,units=timeunits)
                fig_quiver=plotting.quiver_geopot_plot(Udata[2,:,:],Vdata[2,:,:],Phidata[2,:,:]+Phibar, lambdas, mus, timestamp,units=timeunits,minlevel=minlevel,maxlevel=maxlevel) 
                fig_spinup=plotting.spinup_plot(spinupdata, dt,units=timeunits)       
        
        A,B,C,D,E = initial_conditions.ABCDE_init(np.real(newU),np.real(newV),np.real(neweta),np.real(newPhi),mus,I,J)
        
        Adata[2,:,:]=A
        Bdata[2,:,:]=B
        Cdata[2,:,:]=C
        Ddata[2,:,:]=D
        Edata[2,:,:]=E
        
        Amdata[2,:,:]=spectral_transform.fwd_fft_trunc(A, I, M)
        Bmdata[2,:,:]=spectral_transform.fwd_fft_trunc(B, I, M)
        Cmdata[2,:,:]=spectral_transform.fwd_fft_trunc(C, I, M)
        Dmdata[2,:,:]=spectral_transform.fwd_fft_trunc(D, I, M)
        Emdata[2,:,:]=spectral_transform.fwd_fft_trunc(E, I, M)
        
        #rotate all arrays 

        etadata[0:2,:,:]=etadata[1:3,:,:]
        deltadata[0:2,:,:]=deltadata[1:3,:,:]
        Phidata[0:2,:,:]=Phidata[1:3,:,:]

        etamdata[0:2,:,:]=etamdata[1:3,:,:]
        deltamdata[0:2,:,:]=deltamdata[1:3,:,:]
        Phimdata[0:2,:,:]=Phimdata[1:3,:,:]
        
        Udata[0:2,:,:]=Udata[1:3,:,:]
        Vdata[0:2,:,:]=Vdata[1:3,:,:]

        Adata[0:2,:,:]=Adata[1:3,:,:]
        Bdata[0:2,:,:]=Bdata[1:3,:,:]
        Cdata[0:2,:,:]=Cdata[1:3,:,:]
        Ddata[0:2,:,:]=Ddata[1:3,:,:]
        Edata[0:2,:,:]=Edata[1:3,:,:]      
        
        Amdata[0:2,:,:]=Amdata[1:3,:,:]
        Bmdata[0:2,:,:]=Bmdata[1:3,:,:]
        Cmdata[0:2,:,:]=Cmdata[1:3,:,:]
        Dmdata[0:2,:,:]=Dmdata[1:3,:,:]
        Emdata[0:2,:,:]=Emdata[1:3,:,:]   
        
        
        Fdata[0:2,:,:]=Fdata[1:3,:,:]
        Gdata[0:2,:,:]=Gdata[1:3,:,:]
            
        Fmdata[0:2,:,:]=Fmdata[1:3,:,:]
        Gmdata[0:2,:,:]=Gmdata[1:3,:,:]
            
        Phiforcingdata[0:2,:,:]=Phiforcingdata[1:3,:,:]
        Phiforcingmdata[0:2,:,:]= Phiforcingmdata[1:3,:,:]
    
    if verbose==True:
        print('GCM run completed!') 
               