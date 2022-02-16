# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 11:39:39 2021

@author: ek672
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
    

tindex=600
ttoprint=tindex*p.savefreq

dt=50
# etadata=cont.load_and_restore('data\Showman2015\cos_bell_experiment\etadata-k1-0.0002-k2-0.0004',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('data\Showman2015\cos_bell_experiment\deltadata-k1-0.0002-k2-0.0004',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('data\Showman2015\cos_bell_experiment\Phidata-k1-0.0002-k2-0.0004',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0



etadata=cont.load_and_restore('etadata-taudrag-864000-taurad-8640',I)
eta0 = etadata[tindex,:,:]
eta1 = eta0
deltadata = cont.load_and_restore('deltadata-taudrag-864000-taurad-8640',I)
delta0=deltadata[tindex,:,:]
delta1 = delta0
Phidata = cont.load_and_restore('Phidata-taudrag-864000-taurad-8640',I)
Phi0=Phidata[tindex,:,:]
Phi1 = Phi0

# etadata=cont.load_and_restore('data\coriolis\etadata-omega-Earth',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('data\coriolis\deltadata-omega-Earth',I)
# delta0=deltadata[tindex,:,:]
# delta1 = delta0
# Phidata = cont.load_and_restore('data\coriolis\Phidata-omega-Earth',I)
# Phi0=Phidata[tindex,:,:]
# Phi1 = Phi0


# etadata=cont.load_and_restore('data\coriolis\etadata-taudrag-8640000-taurad-8640',I)
# eta0 = etadata[tindex,:,:]
# eta1 = eta0
# deltadata = cont.load_and_restore('data\coriolis\deltadata-taudrag-8640000-taurad-8640',I)
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

etam0=rfl.fwd_fft_trunc(eta0, I, M)
etamn0=rfl.fwd_leg(etam0,J,M,N,Pmn,w)
deltam0=rfl.fwd_fft_trunc(delta0, I, M)
deltamn0=rfl.fwd_leg(deltam0,J,M,N,Pmn,w)

Ucomp,Vcomp=rfl.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
U=np.real(Ucomp)
V=np.real(Vcomp)


spinupdata = np.max(np.sqrt(U[:,:]**2 + V[:,:]**2))
print('RMS winds =' +str(spinupdata)+' m/s')
print('Max zonal winds =' +str(np.max(U[:,:]))+' m/s')
print('Max meridional winds =' +str(np.max(V[:,:]))+' m/s')


plt.plot(lambdas*180/np.pi,U[31,:])
plt.title('Equatorial winds')
plt.show()


#testing_plots.spinup_plot(spinupdata,ttoprint,dt,p.test,p.a1)
# testing_plots.spinup_geopot_plot(Phidata,tmax,dt,test,a1)
testing_plots.zonal_wind_plot(U,mus,ttoprint,50,p.test,p.a1)
#testing_plots.zonal_wind_plot(V,mus,ttoprint,10,p.test,p.a1)
Ubar=np.mean(V,axis=1)
Y = np.arcsin(mus)*180/np.pi #generate latitude

plt.plot(Ubar, Y)

plt.xlabel('mean V, m/s')
plt.ticklabel_format(axis='both', style='sci')
plt.ylabel('latitude')
plt.ticklabel_format(axis='both', style='sci')
plt.title('t='+str(50*p.savefreq*tindex/3600)+' hours, test= Hot Jupiter')
plt.show()
    


testing_plots.quiver_geopot_plot(U,V,Phi0+p.Phibar,lambdas,mus,ttoprint,50,5,p.test,p.a1,np.log10(2*10**6),np.log10(5.5*10**6))
#testing_plots.quiver_temp_plot(U,V,Phi0+p.Phibar,3000,lambdas,mus,ttoprint,200,5,p.test,p.a1,1300,1350)