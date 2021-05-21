# -*- coding: utf-8 -*-
"""
Created on Tue May 11 12:21:28 2021

@author: ek672
"""

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
import forcing
import matplotlib.pyplot as plt
import matplotlib as mpl

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import imageio
import matplotlib.ticker as ticker

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
    


def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)    

def physical_plot(plotdata,mus,lambdas,title):
    """Generates plots for the fourier space with dimensions J, M

    :param plotdata: first variable, JxM+1 
    :type statevar1: float

    
    """
    ## Plot the approximation 
    fig= plt.figure()
    
    # Make data.
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    cp = plt.contourf(X, Y, plotdata,30)
    norm = mpl.colors.Normalize(vmin=0, vmax=5*10**6)
    cb = plt.colorbar(cp,format=ticker.FuncFormatter(fmt),norm=norm)
    #cb.set_clim(0, 5*10**(6))
    plt.title(title)
    plt.show()
    

k1vec=[2*10**(-4),0.002,0.01,0.02]
k2vec=[2*10**(-4),0.002,0.01,0.02]


for i in range(len(k1vec)):
    k1=k1vec[i]
    for j in range(len(k2vec)):
        k2=k2vec[j]
        
        
        Phibar=p.Phibar
        DPhieq=2*p.Phibar #p.DPhieq
        # k1=2*10**(-4) #2 is used in Langton and Laughlin
        # k2=4*10**(-4)
        g=p.g
        pressure=100*250*g/10 #(in Pa)
        R=p.R
        Cp=p.Cp
        sigma=p.sigmaSB
        
        Teq=forcing.DoubleGrayTEqfun(Phibar,DPhieq,lambdas,mus,I,J,k1,k2,pressure,g,R,Cp,sigma)
        Phiic0=R*Teq**0.25
        print(np.max(Teq**0.25))
        # tindex=90
        # ttoprint=tindex*p.savefreq
        
        # etadata=cont.load_and_restore('etadata-k1-0.02-k2-0.002',I)
        # etaic0 = etadata[tindex,:,:]
        # etaic1 = etaic0
        # deltadata = cont.load_and_restore('deltadata-k1-0.0002-k2-0.002',I)
        # deltaic0=deltadata[tindex,:,:]
        # deltaic1 = deltaic0
        # Phidata = cont.load_and_restore('Phidata-k1-0.0002-k2-0.002',I)
        # Phiic0=Phidata[tindex,:,:]
        # Phiic1 = Phiic0
                
        # etam0=rfl.fwd_fft_trunc(etaic0, I, M)
        # print(np.shape(etam0))
        # etamn0=rfl.fwd_leg(etam0,J,M,N,Pmn,w)
        # deltam0=rfl.fwd_fft_trunc(deltaic0, I, M)
        # deltamn0=rfl.fwd_leg(deltam0,J,M,N,Pmn,w)
        
        # Uiccomp,Viccomp=rfl.invrsUV(deltamn0,etamn0,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray)
        Uic=np.zeros((J,I))#np.real(Uiccomp)
        Vic=np.zeros((J,I))
        
        
        
        
        #testing_plots.spinup_plot(spinupdata,tmax,dt,test,a1)
        # testing_plots.spinup_geopot_plot(Phidata,tmax,dt,test,a1)
        # testing_plots.zonal_wind_plot(Uic,mus,ttoprint,10,p.test,p.a1)
        
        # testing_plots.quiver_geopot_plot(Uic,Vic,Phiic0,lambdas,mus,0,10,5,p.test,p.a1,6.1,6.28)
        # testing_plots.quiver_temp_plot(Uic,Vic,Phiic0,p.R,lambdas,mus,0,10,5,p.test,p.a1,500,600)
        # testing_plots.physical_plot(Phiic0, mus, lambdas)
        physical_plot(Teq**0.25, mus, lambdas,'Equilibrium temperature for k1='+str(k1)+', k2='+str(k2))