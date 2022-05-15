# -*- coding: utf-8 -*-
"""
Created on Mon May  2 19:38:05 2022

@author: ek672
"""


import numpy as np


import initial_conditions as ic

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

                
                spinupdata=read_pickle(path, 'spinup-winds')
                
                
                filename='spinup/'+forcingarray[i]+'-Phibar'+phibararray[j]+'-P'+periodarray[l]+'-tau'+tauarray[k]

                write_pickle(filename+'-spinup',spinupdata)


# rms_winds=testing_plots.RMS_winds(p.a, I, J, lambdas, mus, U, V)
# print('RMS winds are '+str(rms_winds))
# plt.show()
# plt.plot(np.arange(len(geopot))*dt/3600,geopot[:,0])
# plt.plot(np.arange(len(geopot))*dt/3600,geopot[:,1])
