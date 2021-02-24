# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 13:18:44 2020

@author: ek672
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:15:17 2020

@author: ek672
"""



import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import math


import initial_conditions as ic
import fft_legendre_trans as rfl
import params as p
import testing_plots
import filters
import tstepping_new
import testing_plots


################ PARAMS #####################
#Spectral parameters
M=p.M #the largest Fourier wave number
N=M #highest degree of the Legendre functions for m=0
K=N+M #highest degree of the Legendre functions
I=p.I #length of array/ number of samples
J=p.J

lambdas=np.linspace(0, 2*np.pi, num=I) #longitudes
[mus,w]=sp.roots_legendre(J) #Gaussian latitudes and weights

Pmn, Hmn= rfl.PmnHmn(J, M, N, mus)


Hmnrecur=np.zeros(np.shape(Hmn))
epsmn=np.zeros((M+2,N+2))

Pmnlarge,Hmnlarge=rfl.PmnHmn(J,M+1,N+1,mus)


for m in range(M+2):
    for n in range(m,N+2):
        epsmn[m,n]=np.sqrt((n**2-m**2)/(4*n**2-1))


for j in range(J):
    for m in range(M+1):
        for n in range(m,N+1):
            Hmnrecur[j,m,n]=-n*epsmn[m,n+1]*Pmnlarge[j,m,n+1]+(n+1)*epsmn[m,n]*Pmnlarge[j,m,n-1]


plt.plot(mus,Hmn[:,1,5])
plt.plot(mus,Hmnrecur[:,1,5])

fmn=np.zeros([M+1,N+1])
fmn[0,1]=p.omega/np.sqrt(0.375)


test,fm=rfl.invrs_leg(fmn,I,J,M,N,Pmn)
f=rfl.invrs_fft(fm,I)
# tstepcoeffmn=tstepping_new.tstepcoeffmn(M,N,p.a)



# #### INITIAL CONDITIONS ####


# etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J)
# Uic,Vic=ic.velocity_init(I,J)
# Aic,Bic,Cic,Dic,Eic=ic.ABCDE_init(Uic,Vic,etaic0,Phiic0,p.mus,I,J)


# #### Get spectral coefficients of eta and delta ####

# etam= rfl.fwd_fft_trunc(etaic0,I,M)
# deltam= rfl.fwd_fft_trunc(deltaic0,I,M)

# etamn=rfl.fwd_leg(etam,J,M,N,Pmn,w)
# deltamn=rfl.fwd_leg(deltam,J,M,N,Pmn,w)

# ##### Compute the U and V from the spectral coeffs of state variables #####

# Unew, Vnew=rfl.invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn)


# ##### Compute the diagnostic relationship eta and delta from U and V #####

# ## from initial conditions ##

# Uicm=rfl.fwd_fft_trunc(Uic,I,M)
# Vicm=rfl.fwd_fft_trunc(Vic,I,M)


# ## initialize arrays for computing the diagnostic relationship
# VargEta=np.zeros(np.shape(Uicm))
# UargEta=np.zeros(np.shape(Uicm))

# UargDelta=np.zeros(np.shape(Uicm))
# VargDelta=np.zeros(np.shape(Uicm))

# for j in range(J):
#     for m in range(M+1):
#         VargEta[j,m]=Vicm[j,m]*(1j)*m/(p.a*(1-mus[j]**2))
#         UargEta[j,m]=Uicm[j,m]/(p.a*(1-mus[j]**2))
        
#         UargDelta[j,m]=Uicm[j,m]*(1j)*m/(p.a*(1-mus[j]**2))
#         VargDelta[j,m]=-Vicm[j,m]/(p.a*(1-mus[j]**2))

# etamn1=rfl.fwd_leg(VargEta,J,M,N,Pmn,w)+rfl.fwd_leg(UargEta,J,M,N,Hmn,w)+fmn
# deltamn1=rfl.fwd_leg(UargDelta,J,M,N,Pmn,w)+rfl.fwd_leg(VargDelta,J,M,N,Hmn,w)

# ### compute physical state variables
# test, etam1=rfl.invrs_leg(etamn1,I,J,M,N,Pmn)
# test, deltam1=rfl.invrs_leg(deltamn1,I,J,M,N,Pmn)

# eta1=rfl.invrs_fft(etam1,I)
# delta1=rfl.invrs_fft(deltam1,I)
 
# ################### from Unew and Vnew #######################
# Unewm=rfl.fwd_fft_trunc(Unew,I,M)
# Vnewm=rfl.fwd_fft_trunc(Vnew,I,M)


# ## initialize arrays for computing the diagnostic relationship
# VargEta=np.zeros(np.shape(Uicm))
# UargEta=np.zeros(np.shape(Uicm))

# UargDelta=np.zeros(np.shape(Uicm))
# VargDelta=np.zeros(np.shape(Uicm))

# for j in range(J):
#     for m in range(M+1):
#         VargEta[j,m]=Vnewm[j,m]*(1j)*m/(p.a*(1-mus[j]**2))
#         UargEta[j,m]=Unewm[j,m]/(p.a*(1-mus[j]**2))
        
#         UargDelta[j,m]=Unewm[j,m]*(1j)*m/(p.a*(1-mus[j]**2))
#         VargDelta[j,m]=-Vnewm[j,m]/(p.a*(1-mus[j]**2))

# etamn2=rfl.fwd_leg(VargEta,J,M,N,Pmn,w)+rfl.fwd_leg(UargEta,J,M,N,Hmn,w)+fmn
# deltamn2=rfl.fwd_leg(UargDelta,J,M,N,Pmn,w)+rfl.fwd_leg(VargDelta,J,M,N,Hmn,w)

# ### compute physical state variables
# test, etam2=rfl.invrs_leg(etamn2,I,J,M,N,Pmn)
# test, deltam2=rfl.invrs_leg(deltamn2,I,J,M,N,Pmn)

# eta2=rfl.invrs_fft(etam2,I)
# delta2=rfl.invrs_fft(deltam2,I)


# ############## STREAM FUNCTION AND VELOCITY POTENTIAL ########################

# ## compute spectral coefficients
# psimn=np.zeros((M+1,N+1))
# ximn=np.zeros((M+1,N+1))
# for m in range(M+1):
#     for n in range(N+1):
#         if n!=0:
#             psimn[m,n]=((p.a)**2/(-n*(n+1)))*(etamn[m,n]-fmn[m,n])
#             ximn[m,n]=((p.a)**2/(-n*(n+1)))*(deltamn[m,n])
    
# test,psim=rfl.invrs_leg(psimn,I,J,M,N,Pmn)
# test,xim=rfl.invrs_leg(ximn,I,J,M,N,Pmn)

# psi=rfl.invrs_fft(psim,I)
# xi=rfl.invrs_fft(psim,I)

# ####### COMPUTE U AND V FROM STREAM FUNCTION INSTEAD #########################

# epsmn=np.zeros((M+1,N+1))
# for m in range(M+1):
#     for n in range(m,N+1):
#         epsmn[m,n]=np.sqrt((n**2-m**2)/(4*n**2-1))
   
# Umn=np.zeros((M+1,N+1))
# Vmn=np.zeros((M+1,N+1))
# for m in range(M+1):
#     for n in range(m,N):
#         Umn[m,n]=(n-1)*epsmn[m,n]*psimn[m,n-1]-(n+2)*epsmn[m,n+1]*psimn[m,n+1]+(1j)*m*ximn[m,n]
#         Vmn[m,n]=-(n-1)*epsmn[m,n]*ximn[m,n-1]+(n+2)*epsmn[m,n+1]*ximn[m,n+1]+(1j)*m*psimn[m,n]
        
# test,Um=rfl.invrs_leg(Umn,I,J,M,N,Pmn)
# test,Vm=rfl.invrs_leg(Vmn,I,J,M,N,Pmn)

# U=rfl.invrs_fft(Um,I)
# V=rfl.invrs_fft(Vm,I)

# ##### PLOTTING ######

# # # U
# # testing_plots.physical_plot(Uic,mus,lambdas)
# # testing_plots.physical_plot(Unew,mus,lambdas)
# # testing_plots.physical_plot(U,mus,lambdas)

# # # eta
# testing_plots.physical_plot(etaic0,mus,lambdas)
# # testing_plots.physical_plot(eta1,mus,lambdas)
# # testing_plots.physical_plot(eta2,mus,lambdas)

# # # delta
# # testing_plots.physical_plot(deltaic0,mus,lambdas)
# # testing_plots.physical_plot(delta1,mus,lambdas)
# # testing_plots.physical_plot(delta2,mus,lambdas)

# testing_plots.physical_plot(psi,mus,lambdas)

