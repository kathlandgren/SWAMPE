#
"""
Spyder Editor

This is the main 2Datmo GCM script (with updates on Legendre and stuff)
"""


## Import statements
# Import python packages
import numpy as np
#import matplotlib.pyplot as plt
#import math

# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import forcing
import filters
import schwartztrauber as S

## Set parameter for type of run
# 0: Debugging, plot at every step
runtype = 0


## Set global parameters
# Spacial and spectral dimensions
# I = p.I 
# J = p.J
M = p.M
# N = p.N
# Length of the run in time steps
tmax = p.tmax
# Sine of the latitude
# mus = p.mus
# Wieghts for integrating
# w = p.w
# lambdas=p.lambdas
g=p.g
taurad=p.taurad
taudrag=p.taudrag
Phibar=p.Phibar
Dheq=p.Dheq
omega=p.omega
a=p.a
a1=p.a1
test=p.test

#normalization for the spherical harmonics
normnum = 1

#colorbar settings for plotting
minlevel=p.minlevel
maxlevel=p.maxlevel


N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
#print(M)
# K4=K4*10**10
#dt=10

forcflag=p.forcflag
diffflag=p.diffflag
zeroflag=p.zeroflag
sigma=filters.sigma(M,N,K4,a,dt)
sigmaPhi=filters.sigmaPhi(M, N, K4, a, dt)
modalflag=p.modalflag
if modalflag==1:
    alpha=p.alpha
# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

SU0, sina, cosa, etaamp, Phiamp =ic.test1_init(a, omega, a1)

## Initialize data arrays 
zetadata=np.zeros((tmax,J,I),dtype=complex)
deltadata=np.zeros((tmax,J,I),dtype=complex)
Phidata=np.zeros((tmax,J,I),dtype=complex)

Udata=np.zeros((tmax,J,I),dtype=complex)
Vdata=np.zeros((tmax,J,I),dtype=complex)

Fdata=np.zeros((tmax,J,I),dtype=complex)
Gdata=np.zeros((tmax,J,I),dtype=complex)
PhiFdata=np.zeros((tmax,J,I),dtype=complex)

spinupdata=np.zeros((tmax,2),dtype=complex)

## Set the initial conditions 

etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,g,omega,a,sina,cosa,etaamp,Phiamp,Phibar,test)
Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)

Phiic0=Phiic0+Phibar
Phiic1=Phiic1+Phibar

f_latlon=S.f_latlon(mus,lambdas,I,J,omega,a1,test)
## Store initial conditions in the data arrays files for easy access
zetadata[0,:,:]=etaic0-f_latlon
zetadata[1,:,:]=etaic1-f_latlon

deltadata[0,:,:]=deltaic0
deltadata[1,:,:]=deltaic1

Phidata[0,:,:]=Phiic0
Phidata[1,:,:]=Phiic1

Udata[0,:,:]=Uic
Vdata[0,:,:]=Vic


#### Forcing ####
heq=forcing.heqfun(Phibar, Dheq, lambdas, mus, I, J,g)
Q=forcing.Qfun(heq, Phiic0, Phibar, taurad, g)
PhiF=Q

if taudrag==-1:
    F=np.divide(np.multiply(-Uic,Q),Phiic0)
    G=np.divide(np.multiply(-Vic,Q),Phiic0)
    
else:
    F=np.divide(np.multiply(-Uic,Q),Phiic0)
    G=np.divide(np.multiply(-Vic,Q),Phiic0)
    F[Q<0]=0
    G[Q<0]=0
    
    F=F-Uic/taudrag
    G=G-Vic/taudrag
    


Fdata[0,:,:]=F
Fdata[1,:,:]=F
Gdata[0,:,:]=G
Gdata[1,:,:]=G
PhiFdata[0,:,:]=PhiF
PhiFdata[1,:,:]=PhiF  
    

# Phiforcingdata[0,:,:]=g*forcing.Qfun(heq, Phiic0, Phibar, taurad, g)
# Phiforcingdata[1,:,:]=g*forcing.Qfun(heq, Phiic1, Phibar, taurad, g)

# Qic=forcing.Qfun(heq, Phiic0, Phibar,taurad,g)

# F=-np.divide(np.multiply(Uic,Qic),(Phiic0+Phibar)/g)
# F[Qic<0]=0
# Fdata[0,:,:]=F
# Fdata[1,:,:]=Fdata[0,:,:]
# G=-np.divide(np.multiply(Vic,Qic),(Phiic0+Phibar)/g)
# G[Qic<0]=0
# Gdata[0,:,:]=G
# Gdata[1,:,:]=Gdata[0,:,:]

# Spin Up calculations
spinupdata[0,0] = np.min(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))
spinupdata[0,1] = np.max(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))


####
# Time stepping
####

## Matrices with coefficients that depend on M, N
nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT=S.mnMATgen(I,J,M,N,mus)

## time -stepping

for t in range(2,tmax):
    print('t='+str(t))
    
    
    delta0=deltadata[t-2,:,:]
    delta1=deltadata[t-1,:,:]
    
    zeta0=zetadata[t-2,:,:]
    zeta1=zetadata[t-1,:,:]

    
    Phi0=Phidata[t-2,:,:]
    Phi1=Phidata[t-1,:,:]
    
    U0=Udata[t-2,:,:]
    V0=Vdata[t-2,:,:]
    
    F0=Fdata[t-2,:,:]
    G0=Gdata[t-2,:,:]
    PhiF0=PhiFdata[t-2,:,:]
        
    
    #forcing 
    # Um=Umdata[t-1,:,:]
    # Vm=Vmdata[t-1,:,:]
    
    # Fm=Fmdata[t-1,:,:]
    # Gm=Gmdata[t-1,:,:]
    

    #PhiFM=Phiforcingmdata[t-1,:,:]    
    
    newdelta, newzeta, newPhi, newU, newV=tstep.tstepping_latlon(test,U0,V0,delta0,delta1,zeta0,zeta1,f_latlon,Phi0,Phi1, w, mus,J,M,nMAT1,nMAT2,nMAT3,mnMAT1,mnMAT2,mnMAT3,mnMAT4,mnMAT5,musMAT,a,dt,Phibar, normnum,forcflag,PhiF0,F0,G0)
    
    #write new data        
    zetadata[t,:,:]=newzeta
    deltadata[t,:,:]=newdelta
    Phidata[t,:,:]=newPhi
    
    if modalflag==1 & t>2:
        print(np.shape(Phidata[t-2:t,:,:]))
        Phidata[t-1,:,:]=filters.modal_splitting(Phidata[t-2:t,:,:],alpha)
        zetadata[t-1,:,:]=filters.modal_splitting(zetadata[t-2:t,:,:],alpha)
        deltadata[t-1,:,:]=filters.modal_splitting(deltadata[t-2:t,:,:],alpha)
        
        
    Udata[t-1,:,:]=newU
    Vdata[t-1,:,:]=newV

    spinupdata[t-1,0] = np.min(np.sqrt(Udata[t-1,:,:]**2 + Vdata[t-1,:,:]**2 ))
    spinupdata[t-1,1] = np.max(np.sqrt(Udata[t-1,:,:]**2 + Vdata[t-1,:,:]**2 ))


    
    #neweta1,newdelta1,etamn1,deltamn1=rfl.diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt)
    
    
    ######## FORCING ############
    # etamdata[t,:,:]=rfl.fwd_fft_trunc(neweta,I,M)
    # deltamdata[t,:,:]=rfl.fwd_fft_trunc(newdelta,I,M)
    # Phimdata[t,:,:]=rfl.fwd_fft_trunc(newPhi,I,M)

    
    Q=forcing.Qfun(heq, newPhi, Phibar,taurad,g)
    PhiF=Q
    
    if taudrag==-1:
        F=np.divide(np.multiply(-newU,Q),Phiic0)
        G=np.divide(np.multiply(-newV,Q),Phiic0)
    
    else:
        F=np.divide(np.multiply(-Uic,Q),Phiic0)
        G=np.divide(np.multiply(-Vic,Q),Phiic0)
        
        F[Q<0]=0
        G[Q<0]=0
        
        F=F-Uic/taudrag
        G=G-Vic/taudrag
    # Phiforcingdata[t,:,:]=g*Q
    # Phiforcingmdata[t,:,:]=rfl.fwd_fft_trunc(Phiforcingdata[t,:,:], I, M)
     
    # F=-np.divide(np.multiply(newU,Q),(newPhi+Phibar)/g)
    # G=-np.divide(np.multiply(newV,Q),(newPhi+Phibar)/g)

    Fdata[t,:,:]=F
    Gdata[t,:,:]=G
    PhiFdata[t,:,:]=PhiF
   
    
    if t%25==0:
        #testing_plots.physical_plot(newPhi,mus,lambdas)
        
        testing_plots.quiver_geopot_plot(newU,newV,newPhi,lambdas,mus,t,dt,6,test,a1,minlevel,maxlevel)
        
        #testing_plots.physical_plot(deltadata[t-1,:,:]-deltadata[t-2,:,:], mus, lambdas)
        #testing_plots.physical_plot(deltadata[t-1,:,:]-deltadata[t-3,:,:], mus, lambdas)
        # testing_plots.physical_plot(neweta-neweta1,mus,lambdas)
        # testing_plots.physical_plot(newV,mus,lambdas)
        # testing_plots.physical_plot(G,mus,lambdas)
        # testing_plots.physical_plot(Q,mus,lambdas)
        testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)
    
   
 
    
    ################# MORE FORCING ##############
    # Fdata[t,:,:]=F
    # Gdata[t,:,:]=G
    
    # Fmdata[t,:,:]=rfl.fwd_fft_trunc(F, I, M)
    # Gmdata[t,:,:]=rfl.fwd_fft_trunc(G, I, M)
    
    
# ####
# #Plotting
# ####
# testing_plots.state_var_compare(etadata[0,:,:], etadata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(deltadata[0,:,:], deltadata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(Phidata[0,:,:]+p.Phibar, Phidata[tmax-1,:,:]+p.Phibar, p.lambdas, p.mus)
# testing_plots.state_var_compare(Udata[0,:,:], Udata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(Vdata[0,:,:], Vdata[tmax-1,:,:], p.lambdas, p.mus)
    
# testing_plots.physical_plot(Phidata[tmax-1,:,:], mus, lambdas)








