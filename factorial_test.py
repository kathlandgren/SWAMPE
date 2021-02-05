# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 18:35:37 2020

@author: ek672
"""

import numpy as np
import math
import testing_plots
import matplotlib.pyplot as plt
import params as p
import scipy.special as sp


M=63
N=63
J=96
scaling_term=np.zeros((M+1,N+1))
for m in range (0,M+1):
    for n in range (m,N+1):
        scaling_term[m,n] = np.sqrt((((2*n)+1)*math.factorial(n-m))/(2*math.factorial(n+m)))
        
testing_plots.spectral_plot(scaling_term)    
plt.plot(scaling_term[0,:])    

def PmnHmn(J,M,N,mus):
    """Calculates the values of associated Legendre polynomials and their 
     derivatives evaluated at Gaussian latitudes (mus) up to wavenumber M 
        
        :param J: number of latitudes
        :type J: int
        :param M: highest wavenumber
        :type M: int
        :param N: highest degree of the Legendre functions for m=0
        :type N: int
        :param mus: Gaussian latitudes
        :type mus: array of float64
        
        :return: Pmn array  (assoc. legendre polynomials), Hmn array (derivatives of Pmn*(1-x^2))
        :rtype: array of float64
        """
    Pmn=np.zeros((J,M+1,N+1))
    Hmn=np.zeros((J,M+1,N+1))
    Pmntemp=np.zeros((J,M+1,N+1))
    Hmntemp=np.zeros((J,M+1,N+1))
    for j in range (0,J):
        Pmntemp[j],Hmntemp[j] = sp.lpmn(M,N,mus[j])
        Hmntemp[j,:,:] = (1-mus[j]**2)*Hmntemp[j,:,:]
        
    for m in range (0,M+1):
        for n in range (m,N+1):
            scaling_term = np.sqrt((((2*n)+1)*math.factorial(n-m))/(2*math.factorial(n+m)))
            Pmn[:,m,n] = scaling_term*Pmntemp[:,m,n]
            Hmn[:,m,n] = scaling_term*Hmntemp[:,m,n]

            if m % 2 == 1:
                if n>0:
                    Pmn[:,m,n]=-Pmn[:,m,n]
                    Hmn[:,m,n]=-Hmn[:,m,n]
                
    return Pmn, Hmn

Pmn, Hmn=PmnHmn(J, M, N, p.mus)
plt.plot(p.w[1]*Pmn[1,0,:])
plt.plot(p.w[1]*Hmn[1,0,:])
# plt.plot(Hmn[0,0,:]/Hmn[0,0,:])
