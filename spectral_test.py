#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 12:45:40 2020

@author: Alice
"""


## Import statements
# Import python packages
import numpy as np
import matplotlib.pyplot as plt
import math

# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots

## Set parameter for type of run
# 0: Debugging, plot at every step
runtype = 0


## Set global parameters
# Spacial and spectral dimensions
I = p.I 
J = p.J
M = 20#p.M
N = 40#p.N
# Length of the run in time steps
tmax = p.tmax
# Sine of the latitude
mus = p.mus
# Wieghts for integrating
w = p.w
# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
#orthogterms=rfl.orthogterms(M, N)

## Check Pmn individually

if runtype==0:
    ## Check to make sure Pmn are the right thing
    # When dividing by the scaling term the even Legendre polynomials (m=0)
    # should have p(-1)=p(1)=1
    ntest = 6
    mtest = 0
    scaling_term = np.sqrt((((2*ntest)+1)*math.factorial(ntest-mtest))/(2*math.factorial(ntest+mtest)))
    plt.plot(mus,Pmn[:,0,ntest]/scaling_term)
    
    ## Check to make sure the Associated Legendre Polynomials are orthonormal
    # First orthogonal check from Wikipedia page on Associated Legendre 
    # polynomials.  In the plot you should get 1 for n=ntest1, m<ntest1 and 0
    # otherwise.
    ntest1 = 15
    orthocheckPmn1 = rfl.fwd_leg(Pmn[:,:,ntest1],J,M,N,Pmn,w)
    testing_plots.spectral_plot(orthocheckPmn1)
    # Second orthogonal check from Wikipedia page on Associated Legendre 
    # polynomials. Note that for m=0, we have set the polynomial to zero so 
    # that we don't get infinity when integrating 1/(1-x^2). In the print out
    # you should get a vector with terms (2*ntest2 +1)/2*m, m>0 and 0 for m=0.
    ntest2 = 20
    Pmnscaled = np.zeros((J, M+1, N+1))
    for j in range (0,J):
        Pmnscaled[j]=Pmn[j,:,:]/(1-mus[j]**2)
    
    Pmnscaled[:,0,:] = 0
    orthocheckPmn2 = rfl.fwd_leg(Pmnscaled[:,:,ntest2],J,M,N,Pmn,w)
    print(orthocheckPmn2[:,ntest2])
    #testing_plots.spectral_plot(orthocheckPmn2)
    
## Check FFT Individually
if runtype==0:
    ## Test the forward Fourier transform. The Fourier transform of a constant 
    # function is delta(m,0) where delta is the Kronecker Delta function.
    constfunc = np.ones((J,I))
    fourier_const = rfl.fwd_fft_trunc(constfunc, I, M)
    testing_plots.fourier_plot(fourier_const,mus)

## Check Forward Legendre Individually
if runtype==0:
    ## Test the forward Legendre transform. 
    testfunc = np.outer(Pmn[:,0,20], np.ones(I))
    testing_plots.physical_plot(testfunc, p.mus, p.lambdas)
    fourier_testfunc = rfl.fwd_fft_trunc(testfunc, I, M)
    testing_plots.fourier_plot(fourier_testfunc,mus)
    legendre_testfunc = rfl.fwd_leg(fourier_testfunc,J,M,N,Pmn,w)
    testing_plots.spectral_plot(legendre_testfunc)

