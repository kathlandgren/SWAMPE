#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 17:15:22 2020

@author: Home
"""


## Import statements
# Import python packages
import numpy as np
#import matplotlib.pyplot as plt
#import math

# Import program packages
import params as p
#import reshapefuns
import initial_conditions as ic
import fft_legendre_trans as rfl
import tstepping_new as tstep
import testing_plots
#import forcing
import filters
import schwartztrauber as S

## Set parameter for type of run
# 0: Debugging, plot at every step
runtype = 0


## Set global parameters
# Spacial and spectral dimensions
# I = p.I 
# J = p.J
M = p.M
# N = p.N
# Length of the run in time steps
tmax = p.tmax
# Sine of the latitude
# mus = p.mus
# Wieghts for integrating
# w = p.w
# lambdas=p.lambdas
g=p.g
taurad=p.taurad
taudrag=p.taudrag
Phibar=p.Phibar
Dheq=p.Dheq
omega=p.omega
a=p.a
a1=p.a1
test=p.test

#normalization for the spherical harmonics
normnum = 1


N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
#print(M)
# K4=K4*10**10
#dt=10

forcflag=p.forcflag
diffflag=p.diffflag
zeroflag=p.zeroflag
sigma=filters.sigma(M,N,K4,a,dt)
sigmaPhi=filters.sigmaPhi(M, N, K4, a, dt)
modalflag=p.modalflag
if modalflag==1:
    alpha=p.alpha
# Associated Legendre Polynomials and their derivatives
Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)

SU0, sina, cosa, etaamp, Phiamp =ic.test1_init(a, omega, a1)

## Set the initial conditions 

etaic0, etaic1, deltaic0, deltaic1, Phiic0, Phiic1=ic.state_var_init(I,J,mus,lambdas,g,omega,a,sina,cosa,etaamp,Phiamp,Phibar,test)
Uic,Vic=ic.velocity_init(I,J,SU0,cosa,sina,mus,lambdas,test)

Phiic0=Phiic0+Phibar
Phiic1=Phiic1+Phibar

f_latlon=S.f_latlon(mus,lambdas,I,J,omega,a1,test)

zetaic0=etaic0-f_latlon
zetaic1=etaic1-f_latlon


## Matrices with coefficients that depend on M, N
nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT=S.mnMATgen(I,J,M,N,mus)


############
# Test 7.2
############

dfromUV, 

zfromUV

Ufromdz,Vfromdz=S.A20_A21(delta1,zeta1,M,nMAT3,mnMAT1,mnMAT2,mnMAT3,w,mus,J,normnum)




