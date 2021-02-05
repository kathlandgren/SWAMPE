#
"""
Spyder Editor

This is the main 2Datmo GCM script (with updates on Legendre and stuff)
"""


## Import statements
# Import python packages
import numpy as np
import matplotlib.pyplot as plt
import math

# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots

## Set parameter for type of run
# 0: Debugging, plot at every step
runtype = 0


## Set global parameters
# Spacial and spectral dimensions
I = p.I 
J = p.J
M = 20#p.M
N = 40#p.N
# Length of the run in time steps
tmax = p.tmax
# Sine of the latitude
mus = p.mus
# Wieghts for integrating
w = p.w
# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
#orthogterms=rfl.orthogterms(M, N)

# if runtype==0:
#     ## Check to make sure Pmn are the right thing
#     # When dividing by the scaling term the even Legendre polynomials (m=0)
#     # should have p(-1)=p(1)=1
#     ntest = 2
#     mtest = 0
#     scaling_term = np.sqrt((((2*ntest)+1)*math.factorial(ntest-mtest))/(2*math.factorial(ntest+mtest)))
#     plt.plot(mus,Pmn[:,0,2]/scaling_term)
    
#     ## Check to make sure the Associated Legendre Polynomials are orthonormal
#     # First orthogonal check from Wikipedia page on Associated Legendre 
#     # polynomials.  In the plot you should get 1 for n=ntest1, m<ntest1 and 0
#     # otherwise.
#     ntest1 = 4
#     orthocheckPmn1 = rfl.fwd_leg(Pmn[:,:,ntest],J,M,N,Pmn,w)
#     testing_plots.spectral_plot(orthocheckPmn1)
#     # Second orthogonal check from Wikipedia page on Associated Legendre 
#     # polynomials. Note that for m=0, we have set the polynomial to zero so 
#     # that we don't get infinity when integrating 1/(1-x^2). In the print out
#     # you should get a vector with terms (2*ntest2 +1)/2*m, m>0 and 0 for m=0.
#     ntest2 = 20
#     Pmnscaled = np.zeros((J, M+1, N+1))
#     for j in range (0,J):
#         Pmnscaled[j]=Pmn[j,:,:]/(1-mus[j]**2)
    
#     Pmnscaled[:,0,:] = 0
#     orthocheckPmn2 = rfl.fwd_leg(Pmnscaled[:,:,20],J,M,N,Pmn,w)
#     print(orthocheckPmn2[:,20])
#     #testing_plots.spectral_plot(orthocheckPmn2)


## Initialize data arrays 
etadata=np.zeros((tmax,J,I),dtype=complex)
deltadata=np.zeros((tmax,J,I),dtype=complex)
Phidata=np.zeros((tmax,J,I),dtype=complex)

etamdata=np.zeros((tmax,J,M+1),dtype=complex)
deltamdata=np.zeros((tmax,J,M+1),dtype=complex)
Phimdata=np.zeros((tmax,J,M+1),dtype=complex)

etamPaddata=np.zeros((tmax,J,I),dtype=complex)
deltamPaddata=np.zeros((tmax,J,I),dtype=complex)
PhimPaddata=np.zeros((tmax,J,I),dtype=complex)

etamndata=np.zeros((tmax,M+1,N+1),dtype=complex)
deltamndata=np.zeros((tmax,M+1,N+1),dtype=complex)
Phimndata=np.zeros((tmax,M+1,N+1),dtype=complex)

Udata=np.zeros((tmax,J,I),dtype=complex)
Vdata=np.zeros((tmax,J,I),dtype=complex)

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

## Set the initial conditions 

etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J)
Uic,Vic=ic.velocity_init(I,J)
Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,p.mus,I,J)



# ####################CLEAN THIS UP!##############################
# #### Get spectral coefficients of eta and delta ####

# etam= rfl.fwd_fft_trunc(etaic0,I,M)
# deltam= rfl.fwd_fft_trunc(deltaic0,I,M)

# etamn=rfl.fwd_leg(etam,J,M,N,Pmn,w)
# deltamn=rfl.fwd_leg(deltam,J,M,N,Pmn,w)

# ##### Compute the U and V from the spectral coeffs of state variables #####

fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fmn[0,1]=p.omega/np.sqrt(0.375)

# tstepcoeffmn=tstep.tstepcoeffmn(M,N,p.a)

# Uic, Vic=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn)

##### Compute the diagnostic relationship eta and delta from U and V #####

## from initial conditions ##

Uicm=rfl.fwd_fft_trunc(Uic,I,M)
Vicm=rfl.fwd_fft_trunc(Vic,I,M)


## initialize arrays for computing the diagnostic relationship
VargEta=np.zeros(np.shape(Uicm))
UargEta=np.zeros(np.shape(Uicm))

UargDelta=np.zeros(np.shape(Uicm))
VargDelta=np.zeros(np.shape(Uicm))

for j in range(J):
    for m in range(M+1):
        VargEta[j,m]=Vicm[j,m]*(1j)*m/(p.a*(1-mus[j]**2))
        UargEta[j,m]=Uicm[j,m]/(p.a*(1-mus[j]**2))
        
        UargDelta[j,m]=Uicm[j,m]*(1j)*m/(p.a*(1-mus[j]**2))
        VargDelta[j,m]=-Vicm[j,m]/(p.a*(1-mus[j]**2))

etamn1=rfl.fwd_leg(VargEta,J,M,N,Pmn,w)+rfl.fwd_leg(UargEta,J,M,N,Hmn,w)+fmn
deltamn1=rfl.fwd_leg(UargDelta,J,M,N,Pmn,w)+rfl.fwd_leg(VargDelta,J,M,N,Hmn,w)

### compute physical state variables
test, etam1=rfl.invrs_leg(etamn1,I,J,M,N,Pmn)
test, deltam1=rfl.invrs_leg(deltamn1,I,J,M,N,Pmn)

etaic0=rfl.invrs_fft(etam1,I)
delta0=rfl.invrs_fft(deltam1,I)

etaic1=etaic0
deltaic1=deltaic0

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

####
#Forward Fourier and Legendre transform
####

## FFT


    

Amdata[0,:,:]=rfl.fwd_fft_trunc(Aic, I, M)
Amdata[1,:,:]=rfl.fwd_fft_trunc(Aic, I, M)

Bmdata[0,:,:]=rfl.fwd_fft_trunc(Bic, I, M)
Bmdata[1,:,:]=rfl.fwd_fft_trunc(Bic, I, M)

Cmdata[0,:,:]=rfl.fwd_fft_trunc(Cic, I, M)
Cmdata[1,:,:]=rfl.fwd_fft_trunc(Cic, I, M)


Dmdata[0,:,:]=rfl.fwd_fft_trunc(Dic, I, M)
Dmdata[1,:,:]=rfl.fwd_fft_trunc(Dic, I, M)

Emdata[0,:,:]=rfl.fwd_fft_trunc(Eic, I, M)
Emdata[1,:,:]=rfl.fwd_fft_trunc(Eic, I, M)

etamdata[0,:,:]=rfl.fwd_fft_trunc(etaic0, I, M)
etamdata[1,:,:]=rfl.fwd_fft_trunc(etaic1, I, M)

deltamdata[0,:,0:M+1]=rfl.fwd_fft_trunc(deltaic0, I, M)
deltamdata[1,:,0:M+1]=rfl.fwd_fft_trunc(deltaic1, I, M)

Phimdata[0,:,0:M+1]=rfl.fwd_fft_trunc(Phiic0, I, M)
Phimdata[1,:,0:M+1]=rfl.fwd_fft_trunc(Phiic1, I, M)

## Forward Legendre

etamndata[0,:,:]=rfl.fwd_leg(etamdata[0,:,:],J,M,N,Pmn,w)
etamndata[1,:,:]=rfl.fwd_leg(etamdata[1,:,:],J,M,N,Pmn,w)

deltamndata[0,:,:]=rfl.fwd_leg(deltamdata[0,:,:],J,M,N,Pmn,w)
deltamndata[1,:,:]=rfl.fwd_leg(deltamdata[1,:,:],J,M,N,Pmn,w)

Phimndata[0,:,:]=rfl.fwd_leg(Phimdata[0,:,:],J,M,N,Pmn,w)
Phimndata[1,:,:]=rfl.fwd_leg(Phimdata[1,:,:],J,M,N,Pmn,w)

####
# Time stepping
####

## time-stepping inputs


fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fmn[0,1]=p.omega/np.sqrt(0.375)

fmn=np.zeros([M+1,N+1])
fmn[0,1]=p.omega/np.sqrt(0.375)

tstepcoeffmn=tstep.tstepcoeffmn(M,N,p.a)
tstepcoeff=tstep.tstepcoeff(J,M)
tstepcoeff2=tstep.tstepcoeff2(J,M)
mJarray=tstep.mJarray(J,M)
narray=tstep.narray(M,N)


## time -stepping

for t in range(2,tmax):
    print('t='+str(t))
    
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
    
    
    newetamn,neweta,newdeltamn,newdelta,newPhimn,newPhi,Unew,Vnew=tstep.tstepping(etam0,etam1,deltam0,deltam1,Phim0,Phim1,I,J,M,N,Am,Bm,Cm,Dm,Em,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,mJarray,narray)

    
    #write new data
    etamndata[t,:,:]=newetamn
    deltamndata[t,:,:]=newdeltamn
    Phimndata[t,:,:]=newPhimn
        
    etadata[t,:,:]=neweta
    deltadata[t,:,:]=newdelta
    Phidata[t,:,:]=newPhi
    
    Udata[t,:,:]=Unew
    Vdata[t,:,:]=Vnew
    
    etamdata[t,:,:]=rfl.fwd_fft_trunc(neweta,I,M)
    deltamdata[t,:,:]=rfl.fwd_fft_trunc(newdelta,I,M)

    
    Phimdata[t,:,:]=rfl.fwd_fft_trunc(newPhi,I,M)
    
    #A,B,C,D,E = ic.ABCDE_init(Unew,Vnew,neweta,newPhi,p.mus,I,J)
    A,B,C,D,E = ic.ABCDE_init(Uic,Vic,neweta,newPhi,p.mus,I,J)
    
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
#testing_plots.state_var_compare(Udata[0,:,:], Udata[tmax-1,:,:], p.lambdas, p.mus)
# testing_plots.state_var_compare(Vdata[0,:,:], Vdata[tmax-1,:,:], p.lambdas, p.mus)
    
testing_plots.physical_plot(Phidata[1,:,:], p.mus, p.lambdas)
testing_plots.physical_plot(Phidata[tmax-1,:,:], p.mus, p.lambdas)

testing_plots.physical_plot(deltadata[0,:,:], p.mus, p.lambdas)
testing_plots.physical_plot(deltadata[tmax-1,:,:], p.mus, p.lambdas)








