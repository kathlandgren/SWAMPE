# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 15:44:34 2021
@author: ezets
"""

import numpy as np


import initial_conditions as ic
import testing_plots
import pickle

M=42
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
#Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
        
forcingarray=['high','medium','low']
phibararray=['1','2','3','4']
tauarray=['0p1','1','10']
periodarray=['1','5','10']

tmax=2660
tmin=2600

assert tmax>tmin, "tmax must be greater than tmin"

trange=np.arange(tmin,tmax)
spacing=10
tspaced=trange[::spacing]



def read_pickle(path,filename):
    infile = open(path+filename,'rb')
    var = pickle.load(infile)
    infile.close()
    return var


def write_pickle(filename,data):
    outfile = open(filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()



num_snapshot=len(tspaced)



paths=np.zeros(108)
for i in range(3):
    for j in range(4):
        for k in range(3):
            for l in range(3):
                path=forcingarray[i]+'_forcing/phibar'+phibararray[j]+'/tau'+tauarray[k]+'period'+periodarray[l]+'/data/'


                #create empty arrays for the variables
                Phidata=np.zeros((num_snapshot,J,I))
                Udata=np.zeros((num_snapshot,J,I))
                Vdata=np.zeros((num_snapshot,J,I))     
                
                for m in range(num_snapshot): 
                    Phidata[m,:,:]=read_pickle(path,'Phi-'+str(tspaced[m]))
                    Udata[m,:,:]=read_pickle(path,'U-'+str(tspaced[m]))
                    Vdata[m,:,:]=read_pickle(path,'V-'+str(tspaced[m]))
    
                # Phi0=np.average(Phidata,axis=0)    
                # U=np.average(Udata,axis=0)  
                # V=np.average(Vdata,axis=0)  
                
                # filename='averages/'+forcingarray[i]+'-Phibar'+phibararray[j]+'-P'+periodarray[l]+'-tau'+tauarray[k]

                # write_pickle(filename+'-Phi',Phi0)
                # write_pickle(filename+'-U',U)
                # write_pickle(filename+'-V',V)
    
                filename='gifs/'+forcingarray[i]+'-Phibar'+phibararray[j]+'-P'+periodarray[l]+'-tau'+tauarray[k]+'.gif'    
                sparseness=4
                frms=5
                #string='taurad10period10-v3.gif'
                test=10
                a1=0
                minlevel=np.min(Phidata)
                maxlevel=np.max(Phidata)   
                dt=36000
                testing_plots.write_quiver_gif(lambdas,mus,Phidata,Udata,Vdata,len(tspaced),frms,filename,sparseness,dt,test,a1,minlevel,maxlevel)
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