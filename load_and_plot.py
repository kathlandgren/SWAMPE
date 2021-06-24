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
    

tindex=291
ttoprint=tindex*p.savefreq

etadata=cont.load_and_restore('data\Showman2015\etadata-k1-0.0002-k2-0.0004',I)
etaic0 = etadata[tindex,:,:]
etaic1 = etaic0
deltadata = cont.load_and_restore('data\Showman2015\deltadata-k1-0.0002-k2-0.0004',I)
deltaic0=deltadata[tindex,:,:]
deltaic1 = deltaic0
Phidata = cont.load_and_restore('data\Showman2015\Phidata-k1-0.0002-k2-0.0004',I)
Phiic0=Phidata[tindex,:,:]
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
testing_plots.zonal_wind_plot(Uic,mus,ttoprint,10,p.test,p.a1)
testing_plots.quiver_geopot_plot(Uic,Vic,Phiic0+p.Phibar,lambdas,mus,ttoprint,10,5,p.test,p.a1,6.1,6.4)
testing_plots.quiver_temp_plot(Uic,Vic,Phiic0+p.Phibar,p.R,lambdas,mus,ttoprint,10,5,p.test,p.a1,700,1400)