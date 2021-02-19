# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 11:59:57 2021

@author: ek672
"""


# Import python packages
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh

# Import program packages
import params as p
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import forcing
import filters
import schwartztrauber as S
import testing_plots




# Set spectral dimensions
M = p.M
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w,normnum=ic.spectral_params(M)


#redefine normnum
normnum=1
# Length of the run in time steps
tmax = p.tmax
#surface gravity
g=p.g
#radiative time scale in Earth days
taurad=p.taurad
#drag time scale in Earth days
taudrag=p.taudrag
#mean geopotential height. In hot Jupiter case, Phibar is the flat nightside thickness
Phibar=p.Phibar
#the difference in radiative-equilibrium thickness between the substellar point and the nightside
DPhieq=p.DPhieq
#rotation rate of the planet, radians per second
omega=p.omega
#planetary radius, meters
a=p.a
#angle for test cases 1 and 2, radians
a1=p.a1
#test case, number
test=p.test

#colorbar settings for plotting
minlevel=p.minlevel
maxlevel=p.maxlevel

#forcing flag
forcflag=p.forcflag

#Coriolis force
f_latlon=ic.f_latlon(mus,lambdas,I,J,omega,a1,test)


## Matrices with coefficients that depend on M, N
nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT=S.mnMATgen(I,J,M,N,mus)

## Set the initial conditions 

#parameters for initializing tests 1 and 2

SU0, sina, cosa, etaamp, Phiamp =ic.test1_init(a, omega, a1)
zetaic0, zetaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,g,omega,Phibar,test,a,sina,cosa,etaamp,Phiamp,f_latlon)
Uic,Vic=ic.velocity_init(I,J,mus,lambdas,test,SU0,cosa,sina)


#generate polynomials
Pmn,Hmn= rfl.PmnHmn(J,M,N,mus)

#initialize c arrays
crmn=np.zeros((J,J))
cimn=np.zeros((J,J))

UmC=np.fft.fft(Uic/I,I,1)[:,0:M+1]
VmC=np.fft.fft(Vic/I,I,1)[:,0:M+1]

UmCos=UmC+np.conj(UmC)
UmSin=(1j)*UmC+np.conj((1j)*UmC)

VmCos=VmC+np.conj(VmC)
VmSin=(1j)*VmC+np.conj((1j)*VmC)

MAT1=np.zeros((J,M+1,N+1))
for m in range(1,M+1):
    for n in range(m,N+1):
        MAT1[:,m,n]=(m/np.sqrt(n*(n+1))) #should it be 1 when 0/0?
        

def fwd_leg(UmCos,UmSin,VmCos,VmSin,J,M,N,Pmn,Hmn):
    cr=np.zeros((J,M+1,N+1),dtype=complex) 
    ci=np.zeros((J,M+1,N+1),dtype=complex)     
    for m in range (0, M+1):
        for j in range (0,J):
            #legterm[j,m,:]=w[j]*(data[j,m])*Pmn[j,m,:] 
            
            cr[j,m,:]=w[j]*(UmCos[j,m]*Hmn[j,m,:]+VmSin[j,m]*np.multiply(MAT1[j,m,:],Pmn[j,m,:])) #scale by (1-mus^2)^(-1)??
            ci[j,m,:]=w[j]*(np.multiply(-MAT1[j,m,:],Pmn[j,m,:])*VmCos[j,m]+UmSin[j,m]*Hmn[j,m,:]) #also scale?
    legcoeffr=np.sum(cr,0)
    legcoeffi=np.sum(ci,0)
    return legcoeffr, legcoeffi


temprmn,tempimn=fwd_leg(UmCos,UmSin,VmCos,VmSin,J,M,N,Pmn,Hmn)

crmn[0:M+1,0:M+1]=temprmn

cimn[0:M+1,0:M+1]=tempimn

zetalm = pysh.expand.SHExpandGLQ(zetaic0, w, mus, norm=normnum,csphase=1,lmax_calc=M)

#zeta=S.A14(temprmn,tempimn,nMAT1,mus,M,J,normnum)

def A14(zetarmn,zetaimn,nMAT1,mus,M,J,normnum):
    
    zetarmnscaled=np.multiply(nMAT1,zetarmn)
    zetaimnscaled=np.multiply(nMAT1,zetaimn)

    #zetamnscaled=np.zeros((2,M+1,M+1))  
    # zetamnscaled=np.zeros((2,J,J))
    # zetamnscaled[0,:M+1,:M+1] = np.transpose(zetarmnscaled)
    # zetamnscaled[1,:M+1,:M+1] = np.transpose(zetaimnscaled)

    #zeta=pysh.expand.MakeGridGLQ(zetamnscaled, mus, lmax=J)   
    datamn=0.5*(zetarmnscaled-(1j)*zetaimnscaled)

    zeta=rfl.invs_sht(datamn,I,J,M,Pmn)

    return zeta

zeta=A14(temprmn,tempimn,nMAT1,mus,M,J,normnum)

#plotting
plt.contourf(lambdas, mus, zetaic0-zeta)
plt.colorbar()
plt.title('zeta error')
plt.show()

#plotting
plt.contourf(lambdas, mus, zetaic0)
plt.colorbar()
plt.title('zeta IC')
plt.show()

plt.contourf(lambdas, mus, zeta)
plt.colorbar()
plt.title('zeta Transform')
plt.show()

