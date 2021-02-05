# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:06:35 2020

@author: ek672
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sp


# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
import filters


import fft_legendre_trans as rfl


M=p.M
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

Pmn=np.zeros((J,M+1,N+1))
Hmn=np.zeros((J,M+1,N+1))
Pmntemp=np.zeros((J,M+1,N+1))
Hmntemp=np.zeros((J,M+1,N+1))
Hmnorig=np.zeros((J,M+1,N+1))

Pmnbig=np.zeros((J,M+2,N+2))
Hmnbig=np.zeros((J,M+2,N+2))
Hmnrecur=np.zeros((J,M+1,N+1))
for j in range (0,J):
    Pmntemp[j],Hmntemp[j] = sp.lpmn(M,N,mus[j])
    Hmnorig[j,:,:]=Hmntemp[j,:,:]
    Hmntemp[j,:,:] = (1-mus[j]**2)*Hmntemp[j,:,:]
    Pmnbig[j],Hmnbig[j] = sp.lpmn(M+1,N+1,mus[j])

def other_scaling(m,n):

    other_scaling=np.sqrt((n**2-m**2)/(4*n**2-1))    
    return other_scaling

for m in range (0,M+1):
    for n in range (m,N+1):


        scaling_term = np.sqrt((((2*n)+1)*math.factorial(n-m))/(2*2*np.pi*math.factorial(n+m)))

        Pmn[:,m,n] = scaling_term*Pmntemp[:,m,n]
        Hmn[:,m,n] = scaling_term*Hmntemp[:,m,n]
        Hmnorig[:,m,n] = scaling_term*Hmnorig[:,m,n]
        #Pmnbig[:,m,n] = scaling_term*Pmnbig[:,m,n]
        
        # if m>0:
        Hmnrecur[:,m,n]=-n*other_scaling(m,n+1)*Pmnbig[:,m,n+1]+(n+1)*other_scaling(m,n)*Pmnbig[:,m,n-1]    

        if m % 2 == 1:
            if n>0:
                Pmn[:,m,n]=-Pmn[:,m,n]
                Hmn[:,m,n]=-Hmn[:,m,n]
                Hmnrecur[:,m,n]=-Hmnrecur[:,m,n]
            
for j in range (0,J):

    Hmnrecur[j,:,:] = (1-mus[j]**2)*Hmnrecur[j,:,:]



testfun=np.zeros((J,I))
truederiv=np.zeros((J,I))

for i in range(I):
    for j in range(J):
        testfun[j,i]=Pmn[j,4,4]#Pmn[j,5,5]+Pmn[j,2,2]
        truederiv[j,i]=Hmnorig[j,4,4] #non-(1-mu**2) Hmn
        #testfun[j,i]=4*np.arcsin(mus[j])#(1-mus[j]**2)**2+4*mus[j]

plt.plot(mus,testfun[:,3])
plt.show()
#transfrom
testm=rfl.fwd_fft_trunc((testfun), I, M)
 #true derivative fourier calculation
tdm=rfl.fwd_fft_trunc((truederiv),I,M) 
testing_plots.fourier_plot(testm, mus) 
testing_plots.fourier_plot(tdm, mus) #true derivative fourier plot

tdmn=rfl.fwd_leg(tdm,J,M,N,Pmn,w) #true derivative spectral 

tdmnHJ=rfl.fwd_leg_HJ(tdm, J, M, N, Pmn)
tdmn2=rfl.fwd_leg_int(tdmnHJ, J, M, N, w)

testing_plots.spectral_plot(tdmn)  #true derivative spectral plot
testing_plots.spectral_plot(tdmn2)
#testderiv
#deriv=rfl.fwd_leg(testm, J, M, N, Hmn,w)
# deriv=-rfl.fwd_leg(testm, J, M, N,Hmnorig,w) #the derivative wrt mu obtained using Hmn --- what we are trying to do
# derivDiff=filters.diffusion(deriv, p.dt, p.a, p.K4, M, N)
# testing_plots.spectral_plot(deriv) #spectral plot of the supposed derivative
# #esting_plots.spectral_plot(derivDiff)
# # print(deriv[0,:])
# test, derivnewm=rfl.invrs_leg(deriv,I,J,M,N,Pmn)
# testing_plots.fourier_plot(derivnewm, mus) #fourier plot of the supposed derivative
# derivnew=-(rfl.invrs_fft(derivnewm, I))

# derivdivide=np.zeros(np.shape(derivnew))

# derivinterim=(derivnew-derivnew[0,:])
# for j in range(J):
#     for i in range(I):

        
#         #derivdivide[j,i]=(derivnew[j,i]/(1-mus[j]**2))
#         derivdivide[j,i]=(derivinterim[j,i]/(1-mus[j]**2)) #works for 0,1

# ##Idea: forward transform the derivative, compare spectral coeffs


# #plt.plot(p.mus,Hmnrecur[:,2,2])
# #plt.show()

# plt.plot(p.mus,testfun[:,3])
# plt.title('testfun')
# plt.show()

# plt.plot(p.mus,truederiv[:,3])
# plt.title('truederiv')
# plt.show()
# # plt.plot(p.mus,derivdivide[:,3])
# # plt.title('derivdivide')
# # plt.show()

# plt.plot(p.mus,derivnew[:,3])
# plt.title('derivnew')
# plt.show()

# plt.plot(p.mus,derivinterim[:,3])
# plt.title('derivinterim')
# plt.show()

# plt.plot(p.mus,Hmn[:,0,1])

# plt.show()

# #plt.plot(p.mus,Hmn[:,5,5]+Hmn[:,2,2])

# nvec=np.arange(N+1)
# coeff=np.multiply(np.multiply(nvec,nvec),np.multiply(nvec+1,nvec+1))-4
# const=p.K4*p.dt/p.a**4
# sigmacoeff1=2*const*coeff
# sigmacoeff=(1+sigmacoeff1)

# #sigmacoeff=(1+2*p.dt*p.K4*coeff/p.a**4)
# sigmas=np.divide(1,sigmacoeff)
