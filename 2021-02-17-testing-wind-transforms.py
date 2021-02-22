
"""
Created on Tue Feb 16 12:16:09 2021

@author: ek672
"""

# Import python packages
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh
import scipy.signal

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
#hyperviscosity filter flag
diffflag=p.diffflag
#hyperviscosity coefficients
sigma=filters.sigma(M,N,K4,a,dt)
sigmaPhi=filters.sigmaPhi(M, N, K4, a, dt)


## Matrices with coefficients that depend on M, N
nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT=S.mnMATgen(I,J,M,N,mus)


G=np.ones((M+1,M+1))

for m in range(0,M+1):
    for n in range(1,M+1):
        G[m,n]=np.arctan(M/6-m/1.5)/np.pi+0.5
        
#Coriolis force
f_latlon=ic.f_latlon(mus,lambdas,I,J,omega,a1,test)

## Set the initial conditions 

#parameters for initializing tests 1 and 2

SU0, sina, cosa, etaamp, Phiamp =ic.test1_init(a, omega, a1)
zetaic0, zetaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,g,omega,Phibar,test,a,sina,cosa,etaamp,Phiamp,f_latlon)
Uic,Vic=ic.velocity_init(I,J,mus,lambdas,test,SU0,cosa,sina)

# Uiccos=np.multiply(Uic,musMAT)
# Viccos=np.multiply(Vic,musMAT)

#testing the forward transfrom
Pmn,Hmn=rfl.PmnHmn(J, M, M, mus)
zetamn=rfl.fwd_sht(zetaic0,I,J,M,Pmn,w)
zetalm = pysh.expand.SHExpandGLQ(zetaic0, w, mus, norm=normnum,csphase=1,lmax_calc=M)

q=2
# zetalm=scipy.signal.decimate(zetalm,q)
# zetacmn0 = np.transpose(zetalm0[0,:-1,:-1])
# zetasmn0 = np.transpose(zetalm0[1,:-1,:-1])

#get wind field from delta, zeta
zetalmPad=np.zeros((2,J,J))
zetalmPad[0,:M+1,:M+1]=zetalm[0,:,:]
zetalmPad[1,:M+1,:M+1]=zetalm[1,:,:]
# zetalmPad[0,:M+1,:32]=zetalm[0,:,:]
# zetalmPad[1,:M+1,:32]=zetalm[1,:,:]


zetanewSH=pysh.expand.MakeGridGLQ(zetalmPad, mus, norm=normnum)


zetanew=rfl.invs_sht(zetamn,I,J,M,Pmn)

# lmax=M+1
# zetacmnminus1 = np.zeros((lmax,lmax))
# zetacmnplus1 = np.zeros((lmax,lmax))
# zetasmnminus1 = np.zeros((lmax,lmax))
# zetasmnplus1 = np.zeros((lmax,lmax))


# zetacmnminus1[:,1:] = np.transpose(zetalm0[0,:-2,:-1])
# zetacmnplus1[:,:] = np.triu(np.transpose(zetalm0[0,1:,:-1]))
# zetasmnminus1[:,1:] = np.transpose(zetalm0[1,:-2,:-1])
# zetasmnplus1[:,:] = np.triu(np.transpose(zetalm0[1,1:,:-1]))


zetatestmn=np.zeros((2,M+1,N+1))
deltatestmn=np.zeros((2,M+1,N+1))

#get delta, zeta coefficient from wind field

brmn,bimn=S.A22_A23(Uic,Vic,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)
 
crmn,cimn=S.A24_A25(Uic,Vic,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum)

# crmn=np.multiply(G,crmn)
# cimn=np.multiply(G,cimn)


crmnnew=np.zeros((M+1,M+1))

for i in range(M+1):
    for j in range(M):
        crmnnew[i,j]=crmn[i,j+1]
        print(crmnnew[i,j])
#crmn=crmnnew

# zetatestmn[0,:-1,:]=crmn[1:,:]
# zetatestmn[1,:-1,:]=cimn[1:,:]

# deltatestmn[0,:-1,:]=brmn[1:,:]
# deltatestmn[1,:-1,:]=bimn[1:,:]
# #inverse transform


#zeta=S.A14(crmn,cimn,nMAT1,mus,M,J,normnum)
zeta=S.A14(crmnnew,cimn,nMAT1,mus,M,J,normnum)
delta=S.A15(brmn,bimn,nMAT1,mus,M,J,normnum)

#get wind field from delta, zeta
Ulm=np.zeros((2,J,J))


#cos 
# Ulm[0,0,0] = mnMAT3[0,0]*zetacmnplus1[0,0]
# Ulm[0,2,0] = -mnMAT3[0,0]*zetacmnplus1[0,0]

U=pysh.expand.MakeGridGLQ(Ulm, mus, norm=normnum)


#U,V=S.A20_A21(deltaic0,zetaic0,M,nMAT3,mnMAT1,mnMAT2,mnMAT3,w,mus,J,normnum)


# #plotting
# plt.contourf(lambdas, mus, Uic-U)
# plt.colorbar()
# plt.title('error')
# plt.show()

# #plotting
# plt.contourf(lambdas, mus, Uic)
# plt.colorbar()
# plt.title('IC')
# plt.show()

# plt.contourf(lambdas, mus, U)
# plt.colorbar()
# plt.title('Transform')
# plt.show()


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

# plt.contourf(lambdas, mus, zetanewSH)
# plt.colorbar()
# plt.title('zeta Transform')
# plt.show()

# #plotting
# plt.contourf(lambdas, mus, deltaic0-delta)
# plt.colorbar()
# plt.title('delta error')
# plt.show()

# #plotting
# plt.contourf(lambdas, mus, deltaic0)
# plt.colorbar()
# plt.title('delta IC')
# plt.show()

# plt.contourf(lambdas, mus, delta)
# plt.colorbar()
# plt.title('delta Transform')
# plt.show()

