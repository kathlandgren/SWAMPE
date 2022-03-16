# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 15:44:34 2021
@author: ezets
"""

import numpy as np
import matplotlib.pyplot as plt

import continuation as cont
import testing_plots
import params as p
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep

import pickle

M=p.M  
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
        
tmin=1900
tmax=2000
trange=np.arange(tmin,tmax)
spacing=10
tspaced=trange[::spacing]

num_snapshot=len(tspaced)
#create empty arrays for the variables
etadata=np.zeros((num_snapshot,J,I))
deltadata=np.zeros((num_snapshot,J,I))
Phidata=np.zeros((num_snapshot,J,I))
Udata=np.zeros((num_snapshot,J,I))
Vdata=np.zeros((num_snapshot,J,I))

for k in range(num_snapshot):

    etadata[k,:,:]=cont.read_pickle('eta-'+str(tspaced[k]))
    deltadata[k,:,:]=cont.read_pickle('delta-'+str(tspaced[k]))   
    Phidata[k,:,:]=cont.read_pickle('Phi-'+str(tspaced[k]))
    Udata[k,:,:]=cont.read_pickle('U-'+str(tspaced[k]))
    Vdata[k,:,:]=cont.read_pickle('V-'+str(tspaced[k]))
    
    
    
sparseness=6
frms=5
string='test-snap.gif'
test=10
a1=0
minlevel=p.minlevel
maxlevel=p.maxlevel    
dt=36000
testing_plots.write_quiver_gif(lambdas,mus,Phidata,Udata,Vdata,10,frms,string,sparseness,dt,test,a1,minlevel,maxlevel)
# tindex=1990
# #etaic0 = cont.load_input('etadata')
# eta0=cont.read_pickle('eta-'+str(tindex))
# eta1 = eta0
# delta0=cont.read_pickle('delta-'+str(tindex))
# delta1 = delta0
# #Phiic0 = cont.load_input('Phidata')
# Phi0=cont.read_pickle('Phi-'+str(tindex))
# Phi1 = Phi0

# U=cont.read_pickle('U-'+str(tindex))

# V=cont.read_pickle('V-'+str(tindex))

# rmswinds=cont.read_pickle('spinup-winds')

# etam0=rfl.fwd_fft_trunc(eta0, I, M)
# etamn0=rfl.fwd_leg(etam0,J,M,N,Pmn,w)
# deltam0=rfl.fwd_fft_trunc(delta0, I, M)
# deltamn0=rfl.fwd_leg(deltam0,J,M,N,Pmn,w)



# # fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
# # fmn[0,1]=p.omega/np.sqrt(0.375)

# f=np.zeros([J,I])
# for i in range(I):
#     for j in range(J):
#         f[j,i]=2*p.omega*mus[j]
# fmn=np.zeros([M+1,N+1]) #TODO make a function in tstep
# fm=rfl.fwd_fft_trunc(f, I, M)
# fmn=rfl.fwd_leg(fm, J, M, N, Pmn, w)

# tstepcoeffmn=tstep.tstepcoeffmn(M,N,p.a)
# tstepcoeff=tstep.tstepcoeff(J,M,dt,mus,p.a)
# tstepcoeff2=tstep.tstepcoeff2(J,M,dt,p.a)
# mJarray=tstep.mJarray(J,M)
# marray=tstep.marray(M, N)
# narray=tstep.narray(M,N)
    

# # Ucomp,Vcomp=rfl.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
# # U=np.real(Ucomp)
# # V=np.real(Vcomp)


# ttoprint=int(tindex*p.savefreq/100)

# dt=10



# # plt.plot(lambdas*180/np.pi,U[31,:])
# # plt.title('Equatorial winds')
# # plt.show()


# #testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)
# # testing_plots.spinup_geopot_plot(Phidata,tmax,dt,test,a1)
# testing_plots.zonal_wind_plot(U,mus,ttoprint,dt,p.test,p.a1)
# #testing_plots.zonal_wind_plot(Upic,mus,ttoprint,200,p.test,p.a1)
# #testing_plots.zonal_wind_plot(V,mus,ttoprint,10,p.test,p.a1)
# # Vbar=np.mean(V,axis=1)
# # Y = np.arcsin(mus)*180/np.pi #generate latitude

# # plt.plot(Vbar, Y)

# # plt.xlabel('mean V, m/s')
# # plt.ticklabel_format(axis='both', style='sci')
# # plt.ylabel('latitude')
# # plt.ticklabel_format(axis='both', style='sci')
# # plt.title('t='+str(200*p.savefreq*tindex/3600)+' hours, test= Hot Jupiter')
# # plt.show()
    

# testing_plots.physical_plot(eta0, mus, lambdas)
# testing_plots.physical_plot(delta0, mus, lambdas)
# testing_plots.quiver_geopot_plot(U,V,Phi0+p.Phibar,lambdas,mus,ttoprint,dt,5,p.test,p.a1,p.minlevel,p.maxlevel)
# #testing_plots.quiver_temp_plot(U,V,Phi0+p.Phibar,3000,lambdas,mus,ttoprint,200,5,p.test,p.a1,1300,1350)

# plt.plot(np.arange(40000)*180/3600,rmswinds[:,1])
# plt.xlabel('time, hours')
# plt.ylable('RMS winds, m/s')
# plt.title('RMS winds for HJ')
# # testing_plots.physical_plot(eta0-etapic, mus, lambdas)
# # testing_plots.physical_plot(delta0-deltapic, mus, lambdas)
# # testing_plots.physical_plot(Phi0-Phipic, mus, lambdas)
# # testing_plots.physical_plot(U-Upic, mus, lambdas)
# # testing_plots.physical_plot(V-Vpic, mus, lambdas)