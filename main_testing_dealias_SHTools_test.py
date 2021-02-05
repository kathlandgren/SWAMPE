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
import pyshtools as pysh

# Import program packages
import params as p
#import reshapefuns
import initial_conditions_SH as ic
import fft_legendre_trans as rfl
import fft_legendre_trans_SHTools as rfl_SH
import tstepping_new as tstep
import testing_plots

## Set parameter for type of run
# 0: Debugging, plot at every step
runtype = 0


## Set global parameters
# Spacial and spectral dimensions
M = p.M
N = M #Need to be the same size for SHTools
lmax_calc = M
I = p.I 
J = p.J

# Length of the run in time steps
tmax = p.tmax
# Omega
omega = p.omega
# a
a = p.a
# dt
dt = p.dt

mus_SH, w_SH = pysh.expand.SHGLQ(N)
# Sine of the latitude
mus = mus_SH# p.mus
# Wieghts for integrating
w = w_SH# p.w
# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl_SH.PmnHmn(J, M, N, mus)
# NOTE: The SHTools Pmn and Spectral space are the transpose of the space we're using
# in our code

#Longitudes
lambdas=np.arange(0, 2*np.pi, 2*np.pi/I) #longitudes


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

etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas)
Uic,Vic=ic.velocity_init(I,J,mus)
Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,mus,I,J)

# testing_plots.physical_plot(Phiic0,p.mus,p.lambdas)


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

# etamdata[0,:,:]=rfl.fwd_fft_trunc(etaic0, I, M)
# etamdata[1,:,:]=rfl.fwd_fft_trunc(etaic1, I, M)

# deltamdata[0,:,0:M+1]=rfl.fwd_fft_trunc(deltaic0, I, M)
# deltamdata[1,:,0:M+1]=rfl.fwd_fft_trunc(deltaic1, I, M)

# Phimdata[0,:,0:M+1]=rfl.fwd_fft_trunc(Phiic0, I, M)
# Phimdata[1,:,0:M+1]=rfl.fwd_fft_trunc(Phiic1, I, M)

## Forward Legendre

# etamndata0 = np.zeros((M+1,N+1),dtype=complex)
# etamndata1 = np.zeros((M+1,N+1),dtype=complex)
# deltamndata0 = np.zeros((M+1,N+1),dtype=complex)
# deltamndata1 = np.zeros((M+1,N+1),dtype=complex)
# Phimndata0 = np.zeros((M+1,N+1),dtype=complex)
# Phimndata1 = np.zeros((M+1,N+1),dtype=complex)

# etamndata0 = pysh.expand.SHExpandGLQC(etadata[0,:,:], w, mus)[0,:,:]
# etamndata0 = etamndata0.transpose()
# etamndata1 = pysh.expand.SHExpandGLQC(etadata[1,:,:], w, mus)[0,:,:]

etamndata[0,:,:] = (pysh.expand.SHExpandGLQC(etadata[0,:,:], w, mus)[0,:,:]).transpose()
etamndata[1,:,:] = (pysh.expand.SHExpandGLQC(etadata[1,:,:], w, mus)[0,:,:]).transpose()

deltamndata[0,:,:] = (pysh.expand.SHExpandGLQC(deltadata[0,:,:], w, mus)[0,:,:]).transpose()
deltamndata[1,:,:] = (pysh.expand.SHExpandGLQC(deltadata[1,:,:], w, mus)[0,:,:]).transpose()

Phimndata[0,:,:] = (pysh.expand.SHExpandGLQC(Phidata[0,:,:], w, mus)[0,:,:]).transpose()
Phimndata[1,:,:] = (pysh.expand.SHExpandGLQC(Phidata[0,:,:], w, mus)[0,:,:]).transpose()

####
# Time stepping
####

## time-stepping inputs

# Transpose
fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fmn[0,1] = omega/np.sqrt(0.375)

tstepcoeffmn=tstep.tstepcoeffmn(M, N, a)
tstepcoeff=tstep.tstepcoeff(J, M, mus, dt, a)
tstepcoeff2=tstep.tstepcoeff2(J,M,dt,a)
mJarray=tstep.mJarray(J,M)
narray=tstep.narray(M,N)


## time -stepping

for t in range(2,tmax):
    print('t='+str(t))
    
    etamn0=etamndata[t-2,:,:]
    etamn1=etamndata[t-1,:,:]
    
    deltamn0=deltamndata[t-2,:,:]
    deltamn1=deltamndata[t-1,:,:]
    
    Phimn0=Phimndata[t-2,:,:]
    Phimn1=Phimndata[t-1,:,:]
    
    Am=Amdata[t-1,:,:]
    Bm=Bmdata[t-1,:,:]
    Cm=Cmdata[t-1,:,:]    
    Dm=Dmdata[t-1,:,:]
    Em=Emdata[t-1,:,:]
    
    Phimn1new,Phimntstep,deltamn1new,deltamntstep,etamn1new,etamntstep = tstep.tstepping_dealias_SH(etamn0,etamn1,deltamn0,deltamn1,Phimn0,Phimn1,I,J,M,N,Am,Bm,Cm,Dm,Em,fmn,Pmn,Hmn,w,tstepcoeff,tstepcoeff2,tstepcoeffmn,mJarray,narray)

    Phimndata[t-1,:,:] = Phimn1new
    etamndata[t-1,:,:] = etamn1new
    deltamndata[t-1,:,:] = deltamn1new
    
    #write new data
    etamndata[t,:,:] = etamntstep[0,:,:].transpose()
    deltamndata[t,:,:] = deltamntstep[0,:,:].transpose()
    Phimndata[t,:,:] = Phimntstep[0,:,:].transpose()
    
    # Create new physical variables
    newPhitstep = pysh.expand.MakeGridGLQC(Phimntstep, mus)
    newetatstep = pysh.expand.MakeGridGLQC(etamntstep, mus)
    newdeltatstep = pysh.expand.MakeGridGLQC(deltamntstep, mus)

    
    Unew,Vnew = rfl_SH.invrsUV(deltamndata[t,:,:],etamndata[t,:,:],fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn)
    # Unew,Vnew=rfl.invrsUV(deltamn1new,etamn1new,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn)
    
    etadata[t,:,:]=newetatstep
    deltadata[t,:,:]=newdeltatstep
    Phidata[t,:,:]=newPhitstep
    
    Udata[t,:,:]=Unew
    Vdata[t,:,:]=Vnew
    
    # etamdata[t,:,:]=rfl.fwd_fft_trunc(newetatstep,I,M)
    # deltamdata[t,:,:]=rfl.fwd_fft_trunc(newdeltatstep,I,M)

    
    # Phimdata[t,:,:]=rfl.fwd_fft_trunc(newPhitstep,I,M)
    
    A,B,C,D,E = ic.ABCDE_init(Unew,Vnew,newetatstep,newPhitstep,mus,I,J)
    
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
    
    









