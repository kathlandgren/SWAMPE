# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 19:12:12 2022

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
        

tmax=2650
tmin=2070
assert tmax>tmin, "tmax must be greater than tmin"
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
    
Phi0=np.average(Phidata,axis=0)    
U=np.average(Udata,axis=0)  
V=np.average(Vdata,axis=0)  

cont.write_pickle('high-Phibar4-P10-tau0p1-Phi',Phi0)
cont.write_pickle('high-Phibar4-P10-tau0p1-U',U)
cont.write_pickle('high-Phibar4-P10-tau0p1-V',V)

rmswinds=cont.read_pickle('spinup-winds')
geopot=cont.read_pickle('spinup-geopot')
    
sparseness=4

# string='taurad10.gif'
test=10
a1=0

dt=30
    
tindex=tmax
#testing_plots.physical_plot(eta0, mus, lambdas)
#testing_plots.physical_plot(delta0, mus, lambdas)
#testing_plots.physical_plot(Phi0, mus, lambdas)
minlevel=(1)
maxlevel=1.5*10**6
# minlevel=np.log10(1.2*10**4)
# maxlevel=np.log10(3.8*10**6)
#testing_plots.quiver_geopot_plot(U,V,Phi0,lambdas,mus,tindex,dt,4,p.test,p.a1,minlevci')
# plt.ylabel('latitude')
# plt.ticklabel_format(axis='both', style='sci')
# plt.title('t='+str(200*p.savefreq*tindex/36el,maxlevel)
#testing_plots.quiver_geopot_plot(U,V,Phi0,lambdas,mus,tindex,dt,4,p.test,p.a1,minlevel,maxlevel)
testing_plots.quiver_geopot_plot(U,V,Phi0,lambdas,mus,tindex,dt,4,p.test,p.a1,np.min(Phi0),np.max(Phi0))
#testing_plots.quiver_temp_plot(U,V,Phi0+p.Phibar,3000,lambdas,mus,ttoprint,200,5,p.test,p.a1,1300,1350)


plt.plot(np.arange(len(rmswinds))*dt/3600,rmswinds[:,1])
plt.xlabel('time, hours')
plt.ylabel('RMS winds, m/s')
plt.title('RMS winds for HJ')
plt.show()

testing_plots.zonal_wind_plot(U,mus,tindex,dt,p.test,p.a1)

rms_winds=testing_plots.RMS_winds(p.a, I, J, lambdas, mus, U, V)
print('RMS winds are '+str(rms_winds))
plt.show()
plt.plot(np.arange(len(geopot))*dt/3600,geopot[:,0])
plt.plot(np.arange(len(geopot))*dt/3600,geopot[:,1])

print('max zonal winds are '+str(np.max(U)))

wind_vec_length=np.sqrt(U**2+V**2)

print('max winds are '+str(np.max(wind_vec_length)))