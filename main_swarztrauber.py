#
"""
Spyder Editor

This is the main SWAMP-E GCM script implementing the Swarztrauber (1996) scheme 7
"""


## Import statements

# Import python packages
import numpy as np
import matplotlib.pyplot as plt

# Import program packages
import params as p
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import forcing
import filters
import schwartztrauber as S


## Set global parameters

# Set spectral dimensions
M = p.M
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w,normnum=ic.spectral_params(M)

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

#Coriolis force
f_latlon=ic.f_latlon(mus,lambdas,I,J,omega,a1,test)

## Set the initial conditions 

#parameters for initializing tests 1 and 2
if test<3:
    SU0, sina, cosa, etaamp, Phiamp =ic.test1_init(a, omega, a1)
    zetaic0, zetaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,g,omega,Phibar,test,a,sina,cosa,etaamp,Phiamp,f_latlon)
    Uic,Vic=ic.velocity_init(I,J,mus,lambdas,test,SU0,cosa,sina)
else:
    zetaic0, zetaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,g,omega,Phibar,test)
    Uic,Vic=ic.velocity_init(I,J,mus,lambdas,test)



## Initialize data arrays 
zetadata=np.zeros((tmax,J,I))
deltadata=np.zeros((tmax,J,I))
Phidata=np.zeros((tmax,J,I))

Udata=np.zeros((tmax,J,I))
Vdata=np.zeros((tmax,J,I))

Fdata=np.zeros((tmax,J,I))
Gdata=np.zeros((tmax,J,I))
PhiFdata=np.zeros((tmax,J,I))

spinupdata=np.zeros((tmax,2))



## Store initial conditions in the data arrays files for easy access
zetadata[0,:,:]=zetaic0
zetadata[1,:,:]=zetaic1

deltadata[0,:,:]=deltaic0
deltadata[1,:,:]=deltaic1

Phidata[0,:,:]=Phiic0
Phidata[1,:,:]=Phiic1

Udata[0,:,:]=Uic
Vdata[0,:,:]=Vic


#### Forcing ####
Phieq=forcing.Phieqfun(Phibar, DPhieq, lambdas, mus, I, J, g)
Q=forcing.Qfun(Phieq, Phiic0, taurad)
#geopotential forcing to be passed to time stepping
PhiF=Q


F,G=forcing.Rfun(Uic, Vic, Q, Phiic0,taudrag)
    

#store forcing initial condition data
Fdata[0,:,:]=F
Fdata[1,:,:]=F
Gdata[0,:,:]=G
Gdata[1,:,:]=G
PhiFdata[0,:,:]=PhiF
PhiFdata[1,:,:]=PhiF  

# Spin Up calculations
spinupdata[0,0] = np.min(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))
spinupdata[0,1] = np.max(np.sqrt(Udata[0,:,:]**2 + Vdata[0,:,:]**2 ))

####
# Time stepping
####

## Matrices with coefficients that depend on M, N
nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT=S.mnMATgen(I,J,M,N,mus)

## time-stepping

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
        
    
    newdelta, newzeta, newPhi, newU, newV=tstep.tstepping_latlon(test,U0,V0,delta0,delta1,zeta0,zeta1,f_latlon,Phi0,Phi1, w, mus,J,M,nMAT1,nMAT2,nMAT3,mnMAT1,mnMAT2,mnMAT3,mnMAT4,mnMAT5,musMAT,a,dt,Phibar, normnum,diffflag,K4,forcflag,PhiF0,F0,G0)
    
    #BAD CODING PRACTICES!!!
    newPhi[0:20,:]=np.mean(newPhi)
    newPhi[76:96,:]=np.mean(newPhi)

    
    
    #write new data        
    #zetadata[t,:,:]=newzeta
    zetadata[t,:,:]=newzeta#zetaic0
    deltadata[t,:,:]=newdelta#deltaic0#newdelta
    Phidata[t,:,:]=newPhi
    


    if modalflag==1:
        if t>2:
            print(t-3)
            print(t-2)
        
            print(np.shape(Phidata[t-2:t+1,:,:]))
            Phidata[t-1,:,:]=filters.modal_splitting(Phidata[t-2:t+1,:,:],alpha)
            zetadata[t-1,:,:]=filters.modal_splitting(zetadata[t-2:t+1,:,:],alpha)
            deltadata[t-1,:,:]=filters.modal_splitting(deltadata[t-2:t+1,:,:],alpha)
    
        
    Udata[t-1,:,:]=newU
    Vdata[t-1,:,:]=newV

    spinupdata[t-1,0] = np.min(np.sqrt(Udata[t-1,:,:]**2 + Vdata[t-1,:,:]**2 ))
    spinupdata[t-1,1] = np.max(np.sqrt(Udata[t-1,:,:]**2 + Vdata[t-1,:,:]**2 ))
    
    ######## FORCING ############
    Q=forcing.Qfun(Phieq, newPhi, taurad)
    #geopotential forcing to be passed to time stepping
    PhiF=Q
    F,G=forcing.Rfun(newU, newV, Q, newPhi,taudrag)
    
    Fdata[t,:,:]=F
    Gdata[t,:,:]=G
    PhiFdata[t,:,:]=PhiF
   
    
    if t%31==0:
        #testing_plots.physical_plot(newPhi,mus,lambdas)
        if test==10:
            PhitoPlot=newPhi-Phibar
        else:
            PhitoPlot=newPhi
            
        testing_plots.quiver_geopot_plot(newU,newV,PhitoPlot,lambdas,mus,t,dt,6,test,a1,minlevel,maxlevel)
        # plt.contourf(lambdas, mus, newzeta)
        # plt.colorbar()
        # plt.title('zeta IC')
        # plt.show()
testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)