# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:01:11 2022

@author: ek672
"""
import numpy as np
import continuation as cont

import plotting

import spectral_transform as st

import initial_conditions as ic
import params as p
# import time_stepping as tstep
# import pickle

M=42
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
Pmn, Hmn = st.PmnHmn(J, M, N, mus)
        

timestamp=1000
#etaic0 = cont.load_input('etadata')
# eta0=cont.read_pickle('eta-'+str(tindex),custompath='C:/Users/ek672/Dropbox/SWAMP-E/SWAMP-E/data/')
# eta1 = eta0
# delta0=cont.read_pickle('delta-'+str(tindex),custompath='C:/Users/ek672/Dropbox/SWAMP-E/SWAMP-E/data/')
# delta1 = delta0
# #Phiic0 = cont.load_input('Phidata')
# Phi0=cont.read_pickle('Phi-'+str(tindex),custompath='C:/Users/ek672/Dropbox/SWAMP-E/SWAMP-E/data/')
# Phi1 = Phi0

# U=cont.read_pickle('U-'+str(tindex),custompath='C:/Users/ek672/Dropbox/SWAMP-E/SWAMP-E/data/')

# V=cont.read_pickle('V-'+str(tindex),custompath='C:/Users/ek672/Dropbox/SWAMP-E/SWAMP-E/data/')

eta, delta, Phi, U, V =cont.load_data(timestamp,'C:/Users/ek672/Dropbox/SWAMP-E/Misc_SWAMP-E_files/data/')

#rmswinds=cont.read_pickle('spinup-winds')
#geopot=cont.read_pickle('spinup-geopot')

etam0=st.fwd_fft_trunc(eta, I, M)
etamn0=st.fwd_leg(etam0,J,M,N,Pmn,w)
deltam0=st.fwd_fft_trunc(delta, I, M)
deltamn0=st.fwd_leg(deltam0,J,M,N,Pmn,w)



# fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
# fmn[0,1]=p.omega/np.sqrt(0.375)

f=np.zeros([J,I])
for i in range(I):
    for j in range(J):
        f[j,i]=2*p.omega*mus[j]
fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
fm=st.fwd_fft_trunc(f, I, M)
fmn=st.fwd_leg(fm, J, M, N, Pmn, w)

# tstepcoeffmn=tstep.tstepcoeffmn(M,N,p.a)
# tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,p.a)
# tstepcoeff2=tstep.tstepcoeff2(J,M,dt,p.a)
# mJarray=tstep.mJarray(J,M)
# marray=tstep.marray(M, N)
# narray=tstep.narray(M,N)
    

# Ucomp,Vcomp=rfl.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
# U=np.real(Ucomp)
# V=np.real(Vcomp)

tmax=p.tmax
ttoprint=int(timestamp*p.savefreq/100)

dt=30#120
print(np.max(V))
# cont.write_pickle('lambdas', lambdas)
# cont.write_pickle('mus', mus)

# plt.plot(lambdas*180/np.pi,U[31,:])
# plt.title('Equatorial winds')
# plt.show()


#testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)
#testing_plots.spinup_geopot_plot(geopot,tindex,dt,p.test,p.a1)
plotting.zonal_wind_plot(U,mus,timestamp,dt,p.test,p.a1)
plotting.zonal_wind_plot(V,mus,timestamp,dt,p.test,p.a1)
#testing_plots.zonal_wind_plot(Upic,mus,ttoprint,200,p.test,p.a1)
#testing_plots.zonal_wind_plot(V,mus,ttoprint,10,p.test,p.a1)
# Vbar=np.mean(V,axis=1)
# Y = np.arcsin(mus)*180/np.pi #generate latitude

# plt.plot(Vbar, Y)

# plt.xlabel('mean V, m/s')
# plt.ticklabel_format(axis='both', style='s00)+' hours, test= Hot Jupiter')
# plt.show()
    

#testing_plots.physical_plot(eta0, mus, lambdas)
#testing_plots.physical_plot(delta0, mus, lambdas)
plotting.physical_plot(Phi, mus, lambdas)
minlevel=(1)
maxlevel=1.5*10**6
# minlevel=np.log10(1.2*10**4)
# maxlevel=np.log10(3.8*10**6)
#testing_plots.quiver_geopot_plot(U,V,Phi0,lambdas,mus,tindex,dt,4,p.test,p.a1,minlevci')
# plt.ylabel('latitude')
# plt.ticklabel_format(axis='both', style='sci')
# plt.title('t='+str(200*p.savefreq*tindex/36el,maxlevel)
plotting.quiver_geopot_plot(U,V,Phi,lambdas,mus,timestamp,dt,4,p.test,p.a1,minlevel,maxlevel)
plotting.quiver_geopot_plot(U,V,Phi,lambdas,mus,timestamp,dt,4,p.test,p.a1,np.min(Phi),np.max(Phi))
#testing_plots.quiver_temp_plot(U,V,Phi0+p.Phibar,3000,lambdas,mus,ttoprint,200,5,p.test,p.a1,1300,1350)



print(np.mean(Phi))