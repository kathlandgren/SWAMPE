# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 11:39:39 2021

@author: ek672
"""
import numpy as np

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
    
    
    


etaic0 = cont.load_input('etadata-k1-0p0002-k2-0p002')
etaic1 = etaic0
deltaic0 = cont.load_input('deltadata-k1-0p0002-k2-0p002')
deltaic1 = deltaic0
Phiic0 = cont.load_input('Phidata-k1-0p0002-k2-0p002')
Phiic1 = Phiic0
        
etam0=rfl.fwd_fft_trunc(etaic0, I, M)
print(np.shape(etam0))
etamn0=rfl.fwd_leg(etam0,J,M,N,Pmn,w)
deltam0=rfl.fwd_fft_trunc(deltaic0, I, M)
deltamn0=rfl.fwd_leg(deltam0,J,M,N,Pmn,w)

Uiccomp,Viccomp=rfl.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
Uic=np.real(Uiccomp)
Vic=np.real(Viccomp)



#testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)
# testing_plots.spinup_geopot_plot(Phidata,tmax,dt,test,a1)
testing_plots.zonal_wind_plot(Uic,mus,p.tmax,10,p.test,p.a1)
testing_plots.quiver_geopot_plot(Uic,Vic,Phiic0+p.Phibar,lambdas,mus,p.tmax,10,5,p.test,p.a1,p.minlevel,p.maxlevel)