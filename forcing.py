# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:28:45 2020

@author: ek672
"""

import numpy as np 

import params as p
import testing_plots
import initial_conditions as ic
import fft_legendre_trans as rfl

def heqfun(Phibar,Dheq,lambdas,mus,I,J,g):
    heqMat=Phibar*np.ones((J,I))#/g #H
    #heqMat=np.zeros((J,I))
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (0,0)
            if -np.pi/2<lambdas[i]<np.pi/2:
                heqMat[j,i]=heqMat[j,i]+Dheq*np.cos(lambdas[i])*(1-mus[j]**2)
            
    return heqMat


def Qfun(heq,Phi,Phibar,taurad,g):
    Q=(1/taurad)*(heq-(Phi))#+Phibar))#/g)
    #Q=(1/taurad)*(heq-(Phi)/g)
    return Q


