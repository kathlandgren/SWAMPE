# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 15:22:08 2021

@author: ek672

This script gets RMS winds from the data
"""

import numpy as np
import matplotlib.pyplot as plt

import continuation as cont
import testing_plots
import params as p
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep

M=p.M  
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
        

fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fmn[0,1]=p.omega/np.sqrt(0.375)

tstepcoeffmn=tstep.tstepcoeffmn(M,N,p.a)
tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,p.a)
tstepcoeff2=tstep.tstepcoeff2(J,M,dt,p.a)
mJarray=tstep.mJarray(J,M)
marray=tstep.marray(M, N)
narray=tstep.narray(M,N)
    

tindex=970
ttoprint=tindex*p.savefreq

dt=200
# etadata=cont.load_and_restore('data\Showman2015\cos_bell_experiment\etadata-k1-0.0002-k2-0.0004',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('data\Showman2015\cos_bell_experiment\deltadata-k1-0.0002-k2-0.0004',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('data\Showman2015\cos_bell_experiment\Phidata-k1-0.0002-k2-0.0004',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0



# etadata=cont.load_and_restore('data\coriolis\etadata-omega-1p6e-05',I)

# deltadata = cont.load_and_restore('data\coriolis\deltadata-omega-1p6e-05',I)


# etadata=cont.load_and_restore('data\coriolis\etadata-omega-Earth',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('data\coriolis\deltadata-omega-Earth',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('data\coriolis\Phidata-omega-Earth',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0


etadata=cont.load_and_restore('etadata-taudrag--1-taurad-864000-dt-200',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
deltadata = cont.load_and_restore('deltadata-taudrag--1-taurad-864000-dt-200',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('data\coriolis\Phidata-taudrag-8640000-taurad-8640',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0
        
# etadata=cont.load_and_restore('Phibar-Phieq-experiments\etadata-Phibar-5000000-DPhieq-7500000',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('Phibar-Phieq-experiments\deltadata-Phibar-5000000-DPhieq-7500000',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('Phibar-Phieq-experiments\Phidata-Phibar-5000000-DPhieq-7500000',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0


# etadata=cont.load_and_restore('data\LiuShowman\etadata-Phibar-4000000-DPhieq-40000.0',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('data\LiuShowman\deltadata-Phibar-4000000-DPhieq-40000.0',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('data\LiuShowman\Phidata-Phibar-4000000-DPhieq-40000.0',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0
rmsdata=np.zeros((2,tindex))
for i in range(tindex):
    delta0=deltadata[i,:,:]
    delta1 = delta0
    eta0 = etadata[i,:,:]
    eta1 = eta0
    
    etam0=rfl.fwd_fft_trunc(eta0, I, M)
    etamn0=rfl.fwd_leg(etam0,J,M,N,Pmn,w)
    deltam0=rfl.fwd_fft_trunc(delta0, I, M)
    deltamn0=rfl.fwd_leg(deltam0,J,M,N,Pmn,w)
    
    Ucomp,Vcomp=rfl.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
    U=np.real(Ucomp)
    V=np.real(Vcomp)
    rmsdata[0,i]=np.min(np.sqrt(U[:,:]**2 + V[:,:]**2))
    rmsdata[1,i]=np.max(np.sqrt(U[:,:]**2 + V[:,:]**2))
    
timeaxis=np.arange(tindex)*dt*p.savefreq/3600

plt.plot(timeaxis,rmsdata[0,:])
plt.plot(timeaxis,rmsdata[1,:])
plt.ylabel('RMS winds, m/s')
plt.xlabel('time, h')