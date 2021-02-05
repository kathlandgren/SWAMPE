
"""
Spyder Editor

This is the main 2Datmo GCM script. 
"""

## import statements
import numpy as np
import scipy.special as sp
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import math
import numba
from numba import jit
import timeit
#from scipy.misc import factorial throws error

#local imports
import params as p
import reshapefuns
import initial_conditions as ic
import rfft_leg as rfl
import tstepping
import testing_plots
import filters

I=p.I #import the parameters
J=p.J
M=p.M
N=p.N
tmax=p.tmax

mus=p.mus
w=p.w
Pmn, Hmn= reshapefuns.PmnHmn(J, M, N, mus)


# initialize data arrays 
etadata=np.zeros((tmax,J,I),dtype=complex)
deltadata=np.zeros((tmax,J,I),dtype=complex)
Phidata=np.zeros((tmax,J,I),dtype=complex)

etamdata=np.zeros((tmax,J,M+1),dtype=complex)
deltamdata=np.zeros((tmax,J,M+1),dtype=complex)
Phimdata=np.zeros((tmax,J,M+1),dtype=complex)

etamPaddata=np.zeros((tmax,J,I),dtype=complex)
deltamPaddata=np.zeros((tmax,J,I),dtype=complex)
PhimPaddata=np.zeros((tmax,J,I),dtype=complex)

etamndata=np.zeros((tmax,M+1,M+1),dtype=complex)
deltamndata=np.zeros((tmax,M+1,M+1),dtype=complex)
Phimndata=np.zeros((tmax,M+1,M+1),dtype=complex)

Udata=np.zeros((tmax,J,I),dtype=complex)
Vdata=np.zeros((tmax,J,I),dtype=complex)

Adata=np.zeros((tmax,J,I),dtype=complex)
Bdata=np.zeros((tmax,J,I),dtype=complex)
Cdata=np.zeros((tmax,J,I),dtype=complex)
Ddata=np.zeros((tmax,J,I),dtype=complex)
Edata=np.zeros((tmax,J,I),dtype=complex)

Amdata=np.zeros((tmax,J,M+1),dtype=complex)
Bmdata=np.zeros((tmax,J,M+1),dtype=complex)
Cmdata=np.zeros((tmax,J,M+1),dtype=complex)
Dmdata=np.zeros((tmax,J,M+1),dtype=complex)
Emdata=np.zeros((tmax,J,M+1),dtype=complex)

etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J)
Uic,Vic=ic.velocity_init(I,J)
Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,p.mustile,I,J)




etadata[0,:,:]=etaic0
etadata[1,:,:]=etaic1

deltadata[0,:,:]=deltaic0
deltadata[1,:,:]=deltaic1

Phidata[0,:,:]=Phiic0
Phidata[1,:,:]=Phiic1

Udata[0,:,:]=Uic
Udata[1,:,:]=Uic

Vdata[0,:,:]=Vic
Vdata[1,:,:]=Vic

Adata[0,:,:]=Aic
Adata[1,:,:]=Aic

Bdata[0,:,:]=Bic
Bdata[1,:,:]=Bic

Cdata[0,:,:]=Cic
Cdata[1,:,:]=Cic

Ddata[0,:,:]=Dic
Ddata[1,:,:]=Dic

Edata[0,:,:]=Eic
Edata[1,:,:]=Eic

####
#Forward Fourier and Legendre transform
####

Amdata[0,:,:]=rfl.fwd_fft_trunc(Aic, I, M)
Bmdata[0,:,:]=rfl.fwd_fft_trunc(Bic, I, M)
Cmdata[0,:,:]=rfl.fwd_fft_trunc(Cic, I, M)
Dmdata[0,:,:]=rfl.fwd_fft_trunc(Dic, I, M)
Emdata[0,:,:]=rfl.fwd_fft_trunc(Eic, I, M)

etamdata[0,:,0:M+1]=rfl.fwd_fft_trunc(etaic0, I, M)
etamdata[1,:,0:M+1]=rfl.fwd_fft_trunc(etaic1, I, M)

deltamdata[0,:,0:M+1]=rfl.fwd_fft_trunc(deltaic0, I, M)
deltamdata[1,:,0:M+1]=rfl.fwd_fft_trunc(deltaic1, I, M)

Phimdata[0,:,0:M+1]=rfl.fwd_fft_trunc(Phiic0, I, M)
Phimdata[1,:,0:M+1]=rfl.fwd_fft_trunc(Phiic1, I, M)


## initialize orthogterms TODO: make a function
orthogterms=np.zeros((M+1,N+1))
for m in range (0,M+1):
    for n in range (m,N+1):
        orthogterms[m,n] = ((2*n) + 1)*math.factorial(n-m)/(2*math.factorial(n+m))


Alegcoeff=rfl.fwd_leg(Amdata[0,:,:],J,M,N,Pmn,w,orthogterms)
Blegcoeff=rfl.fwd_leg(Bmdata[0,:,:],J,M,N,Pmn,w,orthogterms)
Clegcoeff=rfl.fwd_leg(Cmdata[0,:,:],J,M,N,Pmn,w,orthogterms)
Dlegcoeff=rfl.fwd_leg(Dmdata[0,:,:],J,M,N,Pmn,w,orthogterms)
Elegcoeff=rfl.fwd_leg(Emdata[0,:,:],J,M,N,Pmn,w,orthogterms)

etamndata[0,:,:]=rfl.fwd_leg(etamdata[0,:,:],J,M,N,Pmn,w,orthogterms)
etamndata[1,:,:]=rfl.fwd_leg(etamdata[1,:,:],J,M,N,Pmn,w,orthogterms)

deltamndata[0,:,:]=rfl.fwd_leg(deltamdata[0,:,:],J,M,N,Pmn,w,orthogterms)
deltamndata[1,:,:]=rfl.fwd_leg(deltamdata[1,:,:],J,M,N,Pmn,w,orthogterms)

Phimndata[0,:,:]=rfl.fwd_leg(Phimdata[0,:,:],J,M,N,Pmn,w,orthogterms)
Phimndata[1,:,:]=rfl.fwd_leg(Phimdata[1,:,:],J,M,N,Pmn,w,orthogterms)

Phim=rfl.fwd_fft_trunc(Phiic0, I, M)
Phimn=(rfl.fwd_leg(Phim, J, M, N, Pmn, w,orthogterms))
####
#Time-stepping
####

#reshape arrays
AmR=reshapefuns.Xmreshape(Amdata[0,:,:],J,M,N)
BmR=reshapefuns.Xmreshape(Bmdata[0,:,:],J,M,N)
CmR=reshapefuns.Xmreshape(Cmdata[0,:,:],J,M,N)
DmR=reshapefuns.Xmreshape(Dmdata[0,:,:],J,M,N)
EmR=reshapefuns.Xmreshape(Emdata[0,:,:],J,M,N)

#etaic0mR=reshapefuns.Xmreshape(etamdata[0,:,:],J,M,N)
#etaic1mR=reshapefuns.Xmreshape(etamdata[1,:,:],J,M,N)

#deltaic0mR=reshapefuns.Xmreshape(deltamdata[0,:,:],J,M,N)
#deltaic1mR=reshapefuns.Xmreshape(deltamdata[1,:,:],J,M,N)

Phiic0mR=reshapefuns.Xmreshape(Phimdata[0,:,:],J,M,N)
Phiic1mR=reshapefuns.Xmreshape(Phimdata[1,:,:],J,M,N)

wR=reshapefuns.wmureshape(w,J,M,N)
muR=reshapefuns.wmureshape(mus,J,M,N)

MR=reshapefuns.mreshape(J,M,N)
NR=reshapefuns.nreshape(J,M,N)




####
# Time-stepping inputs that don't change
####
fmnIJ=np.zeros([M+1,N+1,I,J])
fmnIJ[0,1,:,:]=p.omega/np.sqrt(0.375)

fmn=np.zeros([M+1,N+1])
fmn[0,1]=p.omega/np.sqrt(0.375)

ms=np.transpose(np.array([np.array([np.array([np.arange(0,M+1),]*(N+1)),]*I),]*J), [3,2,1,0])
ntriu=np.triu(np.array([np.arange(0,N+1),]*(M+1)))
ns=np.transpose(np.array([np.array([ntriu,]*I),]*J), [2,3,1,0])

PmnI=np.transpose(np.array([Pmn,]*I),[2, 3, 0, 1])
HmnI=np.transpose(np.array([Hmn,]*I),[2, 3, 0, 1])

Pmnneg=tstepping.nm(PmnI,ms,ns) #negative m coefficients
Hmnneg=tstepping.nm(HmnI,ms,ns) #negative m coefficients

lambdasMNJ=np.transpose(np.array([np.array([np.array([p.lambdas,]*J),]*(M+1)),]*(N+1)),[1, 0, 3, 2])


tstepcoeffIJ=tstepping.tstepcoeff(M,N,p.a,ns)
tstepcoeffmn=tstepping.tstepcoeffmn(M,N,p.a,ns)





#time-stepping
for t in range(2,tmax):
    print('t='+str(t))
    #reshape arrays into 3-arrays
    eta0mR=reshapefuns.Xmreshape(etamdata[t-2,:,:],J,M,N)
    eta1mR=reshapefuns.Xmreshape(etamdata[t-1,:,:],J,M,N)
    
    delta0mR=reshapefuns.Xmreshape(deltamdata[t-2,:,:],J,M,N)
    delta1mR=reshapefuns.Xmreshape(deltamdata[t-1,:,:],J,M,N)
    
    Phi0mR=reshapefuns.Xmreshape(Phimdata[t-2,:,:],J,M,N)
    Phi1mR=reshapefuns.Xmreshape(Phimdata[t-1,:,:],J,M,N)
    
    AmR=reshapefuns.Xmreshape(Amdata[t-1,:,:],J,M,N)
    BmR=reshapefuns.Xmreshape(Bmdata[t-1,:,:],J,M,N)
    CmR=reshapefuns.Xmreshape(Cmdata[t-1,:,:],J,M,N)
    DmR=reshapefuns.Xmreshape(Dmdata[t-1,:,:],J,M,N)
    EmR=reshapefuns.Xmreshape(Emdata[t-1,:,:],J,M,N)
    
    #etanew, deltanew, Phinew, Unew, Vnew=tstepping.tstep(t, p.dt,p.a,p.K4,M,N,eta0mR,delta0mR,delta1mR,Phi0mR,Phi1mR,Pmn,Hmn,PmnI,HmnI,Pmnneg,Hmnneg,fmnIJ,lambdasMNJ,tstepcoeffIJ,AmR,BmR,CmR,DmR,EmR,MR,NR,muR,wR,ms,ns)
     
    etanew, deltanew, Phinew, Unew, Vnew= tstepping.tstep(t,p.dt,p.a,p.K4,M,N,I,J,eta0mR,delta0mR,delta1mR,Phi0mR,Phi1mR,Pmn,Hmn,fmn,tstepcoeffmn,AmR,BmR,CmR,DmR,EmR,MR,NR,muR,wR)
    
    #write new data
    etamndata[t,:,:]=etanew
    deltamndata[t,:,:]=deltanew
    Phimndata[t,:,:]=Phinew
    
    etamdata[t,:,:],etamPaddata[t,:,:]=rfl.invrs_leg(etanew,I,J,M,N,Pmn)
    deltamdata[t,:,:],deltamPaddata[t,:,:]=rfl.invrs_leg(deltanew,I,J,M,N,Pmn)
    Phimdata[t,:,:],PhimPaddata[t,:,:]=rfl.invrs_leg(Phinew,I,J,M,N,Pmn)
    
    
    etadata[t,:,:]=rfl.invrs_fft(etamPaddata[t,:,:],I);
    deltadata[t,:,:]=rfl.invrs_fft(deltamPaddata[t,:,:],I);
    Phidata[t,:,:]=rfl.invrs_fft(PhimPaddata[t,:,:],I);
    
    Udata[t,:,:]=Unew
    Vdata[t,:,:]=Vnew
    
    # apply filter
    if p.filterflag==1:
        if t>1:
            etadata=filters.filterfun(etadata, t-1, p.alpha)
            deltadata=filters.filterfun(deltadata, t-1, p.alpha)
            Phidata=filters.filterfun(Phidata, t-1, p.alpha)
            
            Udata=filters.filterfun(Udata, t-1, p.alpha)
            Udata=filters.filterfun(Vdata, t-1, p.alpha)
        
    
    Adata[t,:,:]=np.multiply(Unew,etadata[t,:,:])
    Bdata[t,:,:]=np.multiply(Vnew,etadata[t,:,:])
    Cdata[t,:,:]=np.multiply(Unew,Phidata[t,:,:])
    Ddata[t,:,:]=np.multiply(Vnew,Phidata[t,:,:])
    Edata[t,:,:]=np.divide(np.multiply(Unew,Unew)+np.multiply(Vnew,Vnew),2*(1-np.transpose(np.multiply(p.mustile,p.mustile))))
    
    Amdata[t,:,:]=rfl.fwd_fft_trunc(Adata[t,:,:], I, M)
    Bmdata[t,:,:]=rfl.fwd_fft_trunc(Bdata[t,:,:], I, M)
    Cmdata[t,:,:]=rfl.fwd_fft_trunc(Cdata[t,:,:], I, M)
    Dmdata[t,:,:]=rfl.fwd_fft_trunc(Ddata[t,:,:], I, M)
    Emdata[t,:,:]=rfl.fwd_fft_trunc(Edata[t,:,:], I, M)
    
    
####
#Plotting
####
testing_plots.state_var_compare(etadata[0,:,:], etadata[tmax-1,:,:], p.lambdas, p.mus)
testing_plots.state_var_compare(deltadata[0,:,:], deltadata[tmax-1,:,:], p.lambdas, p.mus)
testing_plots.state_var_compare(Phidata[0,:,:], Phidata[tmax-1,:,:], p.lambdas, p.mus)
testing_plots.state_var_compare(Udata[0,:,:], Udata[tmax-1,:,:], p.lambdas, p.mus)
testing_plots.state_var_compare(Vdata[0,:,:], Vdata[tmax-1,:,:], p.lambdas, p.mus)
