#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:18:28 2020

@author: Alice Nadeau
"""

import numpy as np
import scipy.special as sp
import math
import testing_plots
#import pyshtools as pysh


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

def PmnHmnSH(J,M,N,mus):
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
    lmax = M
    temp_size = int((lmax+1)*(lmax+2)/2)
    Pmn=np.zeros((J,M+1,N+1))
    Hmn=np.zeros((J,M+1,N+1))
    Pmntemp=np.zeros((J,temp_size))
    Hmntemp=np.zeros((J,temp_size))
    for j in range (0,J-1):

        Pmntemp[j,:], Hmntemp[j,:] = pysh.legendre.PlmBar_d1(lmax, mus[j])
        Pmntemp[j,:]=0.5*Pmntemp[j,:] #rescale by 1/2 to match our factor
        Hmntemp[j,:] = 0.5*(1-mus[j]**2)*Hmntemp[j,:]
        
    for m in range (0,M+1):
        for n in range (m,N+1):
            Pmn[:,m,n] = Pmntemp[:,pysh.legendre.PlmIndex (n, m)]
            Hmn[:,m,n] = Hmntemp[:,pysh.legendre.PlmIndex (n, m)]
                
    return Pmn, Hmn


def fwd_leg(data,J,M,N,Pmn,w):
    """Calculates the forward legendre transform function

        :param data: input to be transformed (usually output of fft)
        :type data: array of float64 or array of complex128
        :param J: number of latitudes
        :type J: int
        :param M: maximal wave number
        :type M: int
        :param N: highest degree of the Legendre functions for m=0
        :type N: int
        :param Pmn: associated legendre functions evaluated at the Gaussian latitudes mu
        :type Pmn: array of float
        :param w: Gauss Legendre weights
        :type w: array of float64

        :return legcoeff: Legendre coefficients (if data was output from FFT, then legcoeff are the spectral coefficients)
        :rtype legcoeff: array of complex128
        """
    legterm=np.zeros((J,M+1,N+1),dtype=complex) 
    
    for m in range (0, M+1):
        for j in range (0,J):

            legterm[j,m,:]=w[j]*(data[j,m])*Pmn[j,m,:] 
    legcoeff=np.sum(legterm,0)
    return legcoeff

def fwd_leg_HJ(data,J,M,N,PHmn):
    legterm=np.zeros((J,M+1,N+1),dtype=complex) 
    
    for m in range (0, M+1):
        for j in range (0,J):
            legterm[j,m,:]=(data[j,m])*PHmn[j,m,:] 
    return legterm
    
def fwd_leg_int(legterm,J,M,N,w):
    temp=np.zeros((J,M+1,N+1),dtype=complex) 
    
    for m in range (0, M+1):
        for j in range (0,J):
            temp[j,m,:]=w[j]*legterm[j,m,:]
    legcoeff=np.sum(temp,0)
    return legcoeff
    

def fwd_fft_trunc(data,I,M):
    """Calculates and truncates the fast forward fourier transform of the input

        :param data: array of length IxJ (usually the values of variables at lat-long coordinates)
        :type data: array of float64
        :param I: number of longitudes
        :type I: int
        :param M: maximal wave number
        :type M: int

        :return datam: Fourier coefficients 0 through M
        :rtype datam: array of complex128
        """
    datam=np.fft.fft(data/I,I,1)[:,0:M+1]
    
    return datam

def invrs_leg(legcoeff,I,J,M,N,Pmn):
    """Calculates the inverse legendre transfrom function
    
    :param legcoeff: Legendre coefficients
    :type legcoeff: array of complex128
    :param J: number of latitudes
    :type J: int
    :param M: maximal wave number
    :type M: int
    :param Pmn: associated legendre functions and coefficients
    :type Pmn: array of float64
    
    :return: transformed spectral coefficients
    :rtype: array of complex128
    """
    approxXim=np.zeros((J,I),dtype=complex)
    approxXimPos=np.zeros((J,M+1),dtype=complex) 
    approxXimNeg=np.zeros((J,M),dtype=complex) 
    for m in range (0, M+1):

        approxXimPos[:,m]=np.matmul(Pmn[:,m,m:N+1],(legcoeff[m,m:N+1]))
        
        #if m !=M:
        if m !=0:
            negm=-m
            # print(m)
            # print(negm)
            negPmn=((-1)**m)*Pmn[:,m,m:N+1]
            negXileg=((-1)**m)*np.conj(legcoeff[m,m:N+1])
            approxXimNeg[:,negm]=np.matmul(negPmn,negXileg)
    approxXim[:,0:M+1]=approxXimPos

    approxXim[:,I-M:I]=approxXimNeg
    return approxXimPos, approxXim

def invrs_fft(approxXim,I):
    """Calculates the inverse Fourier transform function
    
    :param approxXim: Fourier coefficients
    :type approxXim: array of complex128
    :param I: number of longitudes
    :type I: integer

    :return: long-lat coefficients
    :rtype: array of complex
    """
    approxXinew=np.fft.ifft(I*approxXim,I,1);
    return approxXinew

def invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray):
    
    #do not sum over n=0 according to Hack and Jakob 5.24-5.25
    deltamn[:,0]=0
    etamn[:,0]=0
    # marray=np.zeros((M+1,N+1)) #TODO make this an input, compute once
    # mtemp=np.arange(M+1)
    # for n in range(N+1):
    #     marray[:,n]=mtemp
        
    test,newUm1=invrs_leg((1j)*np.multiply(np.multiply(marray,deltamn),tstepcoeffmn), I,J, M, N, Pmn)
    test,newUm2=invrs_leg(np.multiply(etamn-fmn,tstepcoeffmn), I,J, M, N, Hmn)
    
    test,newVm1=invrs_leg((1j)*np.multiply(np.multiply(marray,etamn-fmn),tstepcoeffmn), I,J, M, N, Pmn)
    test,newVm2=invrs_leg(np.multiply(deltamn,tstepcoeffmn), I,J, M, N, Hmn)
    
    # test,newUm1=invrs_leg((1j)*np.multiply(np.multiply(marray,deltamn),tstepcoeffmn), I,J, M, N, Pmn)
    # test,newUm2=invrs_leg(np.multiply(0,tstepcoeffmn), I,J, M, N, Hmn)

    # test,newVm1=invrs_leg((1j)*np.multiply(np.multiply(marray,0),tstepcoeffmn), I,J, M, N, Pmn)
    # test,newVm2=invrs_leg(np.multiply(deltamn,tstepcoeffmn), I,J, M, N, Hmn)
    
    # test, etam=invrs_leg(etamn, I, J, M, N, Pmn)
    # test, fm=invrs_leg(fmn, I, J, M, N, Pmn)
    # eta=invrs_fft(etam,I)
    # f=invrs_fft(fm,I)
    # import initial_conditions as ic
    # N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(63)
    # testing_plots.physical_plot(eta,mus,lambdas)
    # testing_plots.physical_plot(eta-f,mus,lambdas)
    Unew=-invrs_fft(newUm1-newUm2, I)
    Vnew=-invrs_fft(newVm1+newVm2, I)
    return Unew, Vnew

def diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt):
    coeff=tstepcoeff/(2*dt)
    etacomp1prep=np.multiply(np.multiply(coeff,(1j)*mJarray),Vm)
    etacomp2prep=np.multiply(coeff,Um)
    
    zetamn=fwd_leg(etacomp1prep, J, M, N, Pmn, w)+fwd_leg(etacomp2prep, J, M, N, Hmn, w)
    
    etamn=zetamn+fmn
    
    deltacomp1prep=np.multiply(np.multiply(coeff,(1j)*mJarray),Um)
    deltacomp2prep=np.multiply(coeff,Vm)
    
    deltacomp1=fwd_leg(deltacomp1prep, J, M, N, Pmn, w)
    deltacomp2=fwd_leg(deltacomp2prep, J, M, N, Hmn, w)
    
    deltamn=deltacomp1-deltacomp2
    
    test,newdeltam=invrs_leg(deltamn, I,J, M, N, Pmn)
    newdelta=invrs_fft(newdeltam, I)
    
    test,newetam=invrs_leg(etamn, I,J, M, N, Pmn)
    neweta=invrs_fft(newetam, I)
   
    return neweta,newdelta,etamn,deltamn
    
    