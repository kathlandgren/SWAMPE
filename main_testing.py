#
"""
Spyder Editor

This is the main 2Datmo GCM script (with updates on Legendre and stuff)
"""


## Import statements
# Import python packages
import numpy as np


# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import forcing
import filters

##Set global parameters

# Set spectral dimensions
M = p.M
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

#K4=10**25
dt=150
# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)


# Length of the run in time steps
tmax = p.tmax
#surface gravity
g=p.g
#radiative time scale in Earth days
taurad=p.taurad
#drag time scale in Earth days
taudrag=p.taudrag
#mean geopotential height. In hot Jupiter case, Phibar is the flat nightside thickness
Phibar=p.Phibar
#the difference in radiative-equilibrium thickness between the substellar point and the nightside
DPhieq=p.DPhieq

#rotation rate of the planet, radians per second
omega=p.omega
#planetary radius, meters
a=p.a
#angle for test cases 1 and 2, radians
a1=p.a1
#test case, number
test=p.test

#colorbar settings for plotting
minlevel=p.minlevel
maxlevel=p.maxlevel

#forcing flag
forcflag=p.forcflag
#hyperviscosity filter flag
diffflag=p.diffflag
#hyperviscosity coefficients
sigma=filters.sigma(M,N,K4,a,dt)
sigmaPhi=filters.sigmaPhi(M, N, K4, a, dt)

#flag for anti-aliasing filter as in Hack and Jakob (1992) eq. (4.4)
modalflag=p.modalflag
if modalflag==1:
    alpha=p.alpha
    
## Initialize data arrays 
etadata=np.zeros((tmax,J,I),dtype=complex)
deltadata=np.zeros((tmax,J,I),dtype=complex)
Phidata=np.zeros((tmax,J,I),dtype=complex)

etamdata=np.zeros((tmax,J,M+1),dtype=complex)
deltamdata=np.zeros((tmax,J,M+1),dtype=complex)
Phimdata=np.zeros((tmax,J,M+1),dtype=complex)

etamndata=np.zeros((tmax,M+1,N+1),dtype=complex)
deltamndata=np.zeros((tmax,M+1,N+1),dtype=complex)
Phimndata=np.zeros((tmax,M+1,N+1),dtype=complex)

Udata=np.zeros((tmax,J,I),dtype=complex)
Vdata=np.zeros((tmax,J,I),dtype=complex)

Umdata=np.zeros((tmax,J,M+1),dtype=complex)
Vmdata=np.zeros((tmax,J,M+1),dtype=complex)

Adata=np.zeros((tmax,J,I),dtype=complex)
Bdata=np.zeros((tmax,J,I),dtype=complex)
Cdata=np.zeros((tmax,J,I),dtype=complex)
Ddata=np.zeros((tmax,J,I),dtype=complex)
Edata=np.zeros((tmax,J,I),dtype=complex)

Amdata=np.zeros((tmax,J,M+1),dtype=complex)
Bmdata=np.zeros((tmax,J,M+1),dtype=complex)
Cmdata=np.zeros((tmax,J,M+1),dtype=complex)
Dmdata=np.zeros((tmax,J,M+1),dtype=complex)
Emdata=np.zeros((tmax,J,M+1),dtype=complex)

Fdata=np.zeros((tmax,J,I),dtype=complex)
Gdata=np.zeros((tmax,J,I),dtype=complex)

Fmdata=np.zeros((tmax,J,M+1),dtype=complex)
Gmdata=np.zeros((tmax,J,M+1),dtype=complex)

Phiforcingdata=np.zeros((tmax,J,I),dtype=complex)
Phiforcingmdata=np.zeros((tmax,J,M+1),dtype=complex)

spinupdata=np.zeros((tmax,2),dtype=complex)

## Set the initial conditions 
SU0, sina, cosa, etaamp,Phiamp=ic.test1_init(a, omega, a1)

if test==1:
    etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,test,etaamp,a,sina,cosa,Phibar,Phiamp)
elif test==2:
    etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,test,etaamp,a,sina,cosa,Phibar,Phiamp)   
elif test==10:
    etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,test,etaamp)

Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)
Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)


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



# Spin Up calculations
spinupdata[0,0] = np.min(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))
spinupdata[0,1] = np.max(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))


#### Forcing ####
Phieq=forcing.Phieqfun(Phibar, DPhieq, lambdas, mus, I, J, g)
#Q=forcing.Qfun_with_rampup(Phieq, Phiic0, Phibar,taurad,0,dt)

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


    

Amdata[0,:,:]=rfl.fwd_fft_trunc(Aic, I, M)
Amdata[1,:,:]=Amdata[0,:,:]

Bmdata[0,:,:]=rfl.fwd_fft_trunc(Bic, I, M)
Bmdata[1,:,:]=Bmdata[0,:,:]

Cmdata[0,:,:]=rfl.fwd_fft_trunc(Cic, I, M)
Cmdata[1,:,:]=Cmdata[0,:,:]

Dmdata[0,:,:]=rfl.fwd_fft_trunc(Dic, I, M)
Dmdata[1,:,:]=Dmdata[0,:,:]

Emdata[0,:,:]=rfl.fwd_fft_trunc(Eic, I, M)
Emdata[1,:,:]=Emdata[0,:,:]

etamdata[0,:,:]=rfl.fwd_fft_trunc(etaic0, I, M)
etamdata[1,:,:]=etamdata[0,:,:]

deltamdata[0,:,:]=rfl.fwd_fft_trunc(deltaic0, I, M)
deltamdata[1,:,:]=deltamdata[0,:,:]

Phimdata[0,:,:]=rfl.fwd_fft_trunc(Phiic0, I, M)
Phimdata[1,:,:]=Phimdata[0,:,:]


## Forcing Fourier transform ##

Phiforcingmdata[0,:,:]=rfl.fwd_fft_trunc(Phiforcingdata[0,:,:], I, M)
Phiforcingmdata[1,:,:]=Phiforcingmdata[0,:,:]

Fmdata[0,:,:]=rfl.fwd_fft_trunc(Fdata[0,:,:], I, M)
Fmdata[1,:,:]=Fmdata[0,:,:]

Gmdata[0,:,:]=rfl.fwd_fft_trunc(Gdata[0,:,:], I, M)
Gmdata[1,:,:]=Gmdata[0,:,:]

Umdata[0,:,:]=rfl.fwd_fft_trunc(Udata[0,:,:], I, M)
Umdata[1,:,:]=Umdata[0,:,:]

Vmdata[0,:,:]=rfl.fwd_fft_trunc(Vdata[0,:,:], I, M)
Vmdata[1,:,:]=Vmdata[0,:,:]



## Forward Legendre

etamndata[0,:,:]=rfl.fwd_leg(etamdata[0,:,:],J,M,N,Pmn,w)
etamndata[1,:,:]=etamndata[0,:,:]

deltamndata[0,:,:]=rfl.fwd_leg(deltamdata[0,:,:],J,M,N,Pmn,w)
deltamndata[1,:,:]=deltamndata[0,:,:]

Phimndata[0,:,:]=rfl.fwd_leg(Phimdata[0,:,:],J,M,N,Pmn,w)
Phimndata[1,:,:]=Phimndata[0,:,:]

####
# Time stepping
####

## time-stepping inputs

fmn=ic.coriolismn(M,omega)

tstepcoeffmn=tstep.tstepcoeffmn(M,N,a)
tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,a)
tstepcoeff2=tstep.tstepcoeff2(J,M,dt,a)
mJarray=tstep.mJarray(J,M)
marray=tstep.marray(M, N)
narray=tstep.narray(M,N)


## time -stepping

for t in range(2,tmax):
    # print('t='+str(t))
    
    etam0=etamdata[t-2,:,:]
    etam1=etamdata[t-1,:,:]
    
    deltam0=deltamdata[t-2,:,:]
    deltam1=deltamdata[t-1,:,:]
    
    Phim0=Phimdata[t-2,:,:]
    Phim1=Phimdata[t-1,:,:]
    
    Am=Amdata[t-1,:,:]
    Bm=Bmdata[t-1,:,:]
    Cm=Cmdata[t-1,:,:]    
    Dm=Dmdata[t-1,:,:]
    Em=Emdata[t-1,:,:]
    
    
    #forcing 
    Um=Umdata[t-1,:,:]
    Vm=Vmdata[t-1,:,:]
    
    Fm=Fmdata[t-1,:,:]
    Gm=Gmdata[t-1,:,:]
    

    PhiFM=Phiforcingmdata[t-1,:,:]    
    
    newetamn,neweta,newdeltamn,newdelta,newPhimn,newPhi,newU,newV=tstep.tstepping(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,Fm,Gm,Um,Vm,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,marray,mJarray,narray,PhiFM,dt,a,K4,Phibar,taurad,taudrag,forcflag,diffflag,sigma,sigmaPhi)
    
    
    #write new data
    etamndata[t,:,:]=newetamn
    deltamndata[t,:,:]=newdeltamn
    Phimndata[t,:,:]=newPhimn
        
    etadata[t,:,:]=neweta
    deltadata[t,:,:]=newdelta
    Phidata[t,:,:]=newPhi
    
        
    if modalflag==1:
        if t>2:

            Phidata[t-1,:,:]=filters.modal_splitting(Phidata[t-2:t+1,:,:],alpha)
            etadata[t-1,:,:]=filters.modal_splitting(etadata[t-2:t+1,:,:],alpha)
            deltadata[t-1,:,:]=filters.modal_splitting(deltadata[t-2:t+1,:,:],alpha)
    

    if test==1:
        newU=Uic
        newV=Vic


    
    Udata[t,:,:]=newU
    Vdata[t,:,:]=newV
    
    spinupdata[t-1,0] = np.min(np.sqrt(Udata[t-1,:,:]**2 + Vdata[t-1,:,:]**2 ))
    spinupdata[t-1,1] = np.max(np.sqrt(Udata[t-1,:,:]**2 + Vdata[t-1,:,:]**2 ))
    
    if spinupdata[t-1,1]>8000:
        print('Time stepping stopped due to wind blow up. Max RMS winds = '+str(spinupdata[t-1,1]))
        t=tmax

    
    Umdata[t,:,:]=rfl.fwd_fft_trunc(newU,I,M)
    Vmdata[t,:,:]=rfl.fwd_fft_trunc(newV,I, M)
    
    Um=Umdata[t,:,:]
    Vm=Vmdata[t,:,:]
    
    neweta1,newdelta1,etamn1,deltamn1=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)
    #print('Diagnostic eta - timestepping eta '+str(np.max(neweta1-neweta)))
    
    etamdata[t,:,:]=rfl.fwd_fft_trunc(neweta,I,M)
    deltamdata[t,:,:]=rfl.fwd_fft_trunc(newdelta,I,M)
    Phimdata[t,:,:]=rfl.fwd_fft_trunc(newPhi,I,M)
    
        
    ######## FORCING ############
    # Q=forcing.Qfun_with_rampup(Phieq, newPhi,Phibar, taurad,t,dt)
    Q=forcing.Qfun(Phieq, newPhi,Phibar, taurad)
    #geopotential forcing to be passed to time stepping
    PhiF=Q
    F,G=forcing.Rfun(newU, newV, Q, newPhi,Phibar,taudrag)
    
    Fdata[t,:,:]=F
    Gdata[t,:,:]=G
    
    Fmdata[t,:,:]=rfl.fwd_fft_trunc(F, I, M)
    Gmdata[t,:,:]=rfl.fwd_fft_trunc(G, I, M)
    
    Phiforcingdata[t,:,:]=PhiF
    Phiforcingmdata[t,:,:]=rfl.fwd_fft_trunc(Phiforcingdata[t,:,:], I, M)  
    
    if t%5==0:
        print('t='+str(t))
        #testing_plots.physical_plot(newPhi,mus,lambdas)
        
        #testing_plots.quiver_geopot_plot(newU,newV,newPhi,lambdas,mus,t,6)
        # testing_plots.physical_plot(newdelta-newdelta1, mus, lambdas)
        # testing_plots.physical_plot(neweta-neweta1,mus,lambdas)
        # testing_plots.physical_plot(newV,mus,lambdas)
        # testing_plots.physical_plot(G,mus,lambdas)
        # testing_plots.physical_plot(Q,mus,lambdas)
        
        testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)
        
        testing_plots.quiver_geopot_plot(Udata[t,:,:],Vdata[t,:,:],Phidata[t,:,:]+Phibar,lambdas,mus,t,dt,6,test,a1,minlevel,maxlevel)
        
        # plt.contourf(lambdas, mus, newzeta)
        # plt.colorbar()
        # plt.title('zeta IC')
        # plt.show()

    
    A,B,C,D,E = ic.ABCDE_init(newU,newV,neweta,newPhi,mus,I,J)
    
    Adata[t,:,:]=A
    Bdata[t,:,:]=B
    Cdata[t,:,:]=C
    Ddata[t,:,:]=D
    Edata[t,:,:]=E
    
    Amdata[t,:,:]=rfl.fwd_fft_trunc(A, I, M)
    Bmdata[t,:,:]=rfl.fwd_fft_trunc(B, I, M)
    Cmdata[t,:,:]=rfl.fwd_fft_trunc(C, I, M)
    Dmdata[t,:,:]=rfl.fwd_fft_trunc(D, I, M)
    Emdata[t,:,:]=rfl.fwd_fft_trunc(E, I, M)
    
    

    
# ####
# #Plotting
# ####
# testing_plots.state_var_compare(etadata[0,:,:], etadata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(deltadata[0,:,:], deltadata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(Phidata[0,:,:]+p.Phibar, Phidata[tmax-1,:,:]+p.Phibar, p.lambdas, p.mus)
# testing_plots.state_var_compare(Udata[0,:,:], Udata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(Vdata[0,:,:], Vdata[tmax-1,:,:], p.lambdas, p.mus)
    
testing_plots.physical_plot(Phidata[tmax-1,:,:], mus, lambdas)








