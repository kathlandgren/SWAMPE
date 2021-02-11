#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:08:14 2020

@author: Home
"""


import numpy as np
#import scipy.special as sp
#import math
#import testing_plots
import pyshtools as pysh

def mnMATgen(I,J,M,N,mus):
    
    nMAT1=np.zeros((M+1,N+1))
    nMAT2=np.zeros((M+1,N+1))
    nMAT3=np.zeros((M+2,N+2))
    mnMAT1=np.zeros((M+1,N+1))
    mnMAT2=np.zeros((M+1,N+1))
    mnMAT3=np.zeros((M+1,N+1))
    mnMAT4=np.zeros((M+1,N+1))
    mnMAT5=np.zeros((M+1,N+1))
    musMAT=np.zeros((J,I))
    
    nrange=np.arange(N+1)
    nplus1=nrange+1
    
    nnplus1=np.multiply(nrange, nplus1)
    nnplus1sqrt=np.sqrt(np.multiply(nrange, nplus1))

    for m in range(M+1):
        nMAT1[m,:]=nnplus1sqrt
        
    for m in range(M+1):
        nMAT2[m,:]=nnplus1
        
    for m in range(1,M+1):
        nrangeBIG=np.arange(N+2)
        nplus1BIG=nrangeBIG+1
        nnplus1sqrtBIG=np.sqrt(np.multiply(nrangeBIG, nplus1BIG))
        nMAT3[m,1:]=1/nnplus1sqrtBIG[1:]
        
    for m in range(1,M+1):
        mnMAT1[m,1:]=m/nnplus1sqrt[1:]
        
    for m in range(M+1):
        for n in range(2,N+1):
            mnMAT2[m,n]=(n-1)*(n-m)/((2*(n-1)+1)*np.sqrt((n-1)*n))
    mnMAT2=np.triu(mnMAT2)
    
    for m in range(M+1):
        for n in range(N+1):
            mnMAT3=(n+2)*(n+1+m)/(np.sqrt((n+1)*(n+2))*(2*(n+1)+1))
            
    for m in range(1,M+1):
        for n in range(1,N+1):
            mnMAT4=(n+1)*(n+m)/(np.sqrt(n*(n+1))*(2*n+1))
    
    for m in range(1,M+1):
        for n in range(1,N+1):
            mnMAT5=(n*(n-m+1))/(np.sqrt(n*(n+1))*(2*n+1))
    #mnMAT5=np.triu(mnMAT5)        
    
    for i in range(I):
        musMAT[:,i]=np.divide(1,(1-mus**2))
    
    return nMAT1, nMAT2, nMAT3, mnMAT1, mnMAT2, mnMAT3, mnMAT4, mnMAT5, musMAT

def A14(zetarmn,zetaimn,nMAT1,mus,M,J,normnum):
    
    zetarmnscaled=np.multiply(nMAT1,zetarmn)
    zetaimnscaled=np.multiply(nMAT1,zetaimn)

    #zetamnscaled=np.zeros((2,M+1,M+1))  
    zetamnscaled=np.zeros((2,J,J))
    zetamnscaled[0,:M+1,:M+1] = np.transpose(zetarmnscaled)
    zetamnscaled[1,:M+1,:M+1] = np.transpose(zetaimnscaled)

    #zeta=pysh.expand.MakeGridGLQ(zetamnscaled, mus, lmax=J)   
    zeta=pysh.expand.MakeGridGLQ(zetamnscaled, mus, norm=normnum)
    
    return zeta


def A15(deltarmn,deltaimn,nMAT1,mus,M,J,normnum):
    
    deltarmnscaled=np.multiply(nMAT1,deltarmn)
    deltaimnscaled=np.multiply(nMAT1,deltaimn)

    #deltamnscaled=np.zeros((2,M+1,M+1))  
    deltamnscaled=np.zeros((2,J,J))  
    deltamnscaled[0,:M+1,:M+1] = np.transpose(deltarmnscaled)
    deltamnscaled[1,:M+1,:M+1] = np.transpose(deltaimnscaled)
    
    delta=-pysh.expand.MakeGridGLQ(deltamnscaled, mus, norm=normnum)
    
    return delta

def A20_A21(delta,zeta,M,nMAT3,mnMAT1,mnMAT2,mnMAT3,w,mus,J,normnum):
    """Calculates equations A20 and A21 in Schwartztrauber 1996
        
        :param M: highest wavenumber
        :type M: int
        :param mus: Gaussian latitudes
        :type mus: array of float64
        
        :return: bcmn
        :rtype: array of float32
        :return: bsmn
        :rtype: array of float32
        """
    lmax = M+1
    
        
    deltacmnminus1 = np.zeros((lmax,lmax))
    deltacmnplus1 = np.zeros((lmax,lmax))
    deltasmnminus1 = np.zeros((lmax,lmax))
    deltasmnplus1 = np.zeros((lmax,lmax))
    
    zetacmnminus1 = np.zeros((lmax,lmax))
    zetacmnplus1 = np.zeros((lmax,lmax))
    zetasmnminus1 = np.zeros((lmax,lmax))
    zetasmnplus1 = np.zeros((lmax,lmax))
    
    Umn=np.zeros((2,lmax,lmax))
    Vmn=np.zeros((2,lmax,lmax))
    
    Ulm=np.zeros((2,J,J))
    Vlm=np.zeros((2,J,J))
    #temp_size = int((lmax+1)*(lmax+2)/2)
    
    # Ycmnminus1 = np.zeros((lmax+1,lmax+1))
    # Ycmnplus1 = np.zeros((lmax+1,lmax+1))
    # Ysmnminus1 = np.zeros((lmax+1,lmax+1))
    # Ysmnplus1 = np.zeros((lmax+1,lmax+1))
    # brmn = np.zeros((lmax+1,lmax+1))
    # bimn = np.zeros((lmax+1,lmax+1))
    
    #7.1 in Swarztrauber
    deltalm = -np.multiply(nMAT3,pysh.expand.SHExpandGLQ(delta, w, mus,norm=normnum,csphase=1,lmax_calc=lmax))
    deltacmn = np.transpose(deltalm[0,:-1,:-1]) #cosine coeff
    deltasmn = np.transpose(deltalm[1,:-1,:-1]) #sine coeff
    
    #7.1 in Swarztrauber
    zetalm = np.multiply(nMAT3,pysh.expand.SHExpandGLQ(zeta, w, mus, norm=normnum,csphase=1,lmax_calc=lmax))
    zetacmn = np.transpose(zetalm[0,:-1,:-1])
    zetasmn = np.transpose(zetalm[1,:-1,:-1])


    deltacmnminus1[:,1:] = np.transpose(deltalm[0,:-2,:-1])
    deltacmnplus1[:,:] = np.triu(np.transpose(deltalm[0,1:,:-1]))
    deltasmnminus1[:,1:] = np.transpose(deltalm[1,:-2,:-1])
    deltasmnplus1[:,:] = np.triu(np.transpose(deltalm[1,1:,:-1]))
    
    zetacmnminus1[:,1:] = np.transpose(zetalm[0,:-2,:-1])
    zetacmnplus1[:,:] = np.triu(np.transpose(zetalm[0,1:,:-1]))
    zetasmnminus1[:,1:] = np.transpose(zetalm[1,:-2,:-1])
    zetasmnplus1[:,:] = np.triu(np.transpose(zetalm[1,1:,:-1]))
    
    
    # Ucmnminus1[1:,:] = np.transpose(Ylm[0,:-1,:])
    # Ucmnplus1[:-1,:] = np.transpose(Ylm[0,1:,:])
    # smnminus1[1:,:] = np.transpose(Ylm[1,:-1,:])
    # Ysmnplus1[:-1,:] = np.transpose(Ylm[1,1:,:])
    
    #cos
    Umn[0,:,:] = np.multiply(mnMAT1,deltasmn) - np.multiply(mnMAT2,zetacmnminus1) + np.multiply(mnMAT3,zetacmnplus1)
    #sin
    Umn[1,:,:] = -np.multiply(mnMAT1,deltacmn) - np.multiply(mnMAT2,zetasmnminus1) + np.multiply(mnMAT3,zetasmnplus1)

    #cos
    Vmn[0,:,:] = -np.multiply(mnMAT1,zetasmn) - np.multiply(mnMAT2,deltacmnminus1) + np.multiply(mnMAT3,deltacmnplus1)
    #sin
    Vmn[1,:,:] = np.multiply(mnMAT1,zetacmn) - np.multiply(mnMAT2,deltasmnminus1) + np.multiply(mnMAT3,deltasmnplus1)
    
    Ulm[0,:lmax,:lmax]=np.transpose(Umn[0,:,:])
    Ulm[1,:lmax,:lmax]=np.transpose(Umn[1,:,:])
    
    
    Vlm[0,:lmax,:lmax]=np.transpose(Vmn[0,:,:])
    Vlm[1,:lmax,:lmax]=np.transpose(Vmn[1,:,:])
    ## try lmax option in makegrid, if it fails, let's pad it as we did before to satisfy the 2/3 rule
    
    # print(np.shape(Ulm))
    # print(np.shape(mus))
    # Upad=np.zeros((2,J,J))
    # Upad[:,0:64,0:64]=Ulm
    
    # Vpad=np.zeros((2,J,J))
    # Vpad[:,0:,0:64]=Ulm
    # print(np.shape(Utest))
    U=pysh.expand.MakeGridGLQ(Ulm, mus, norm=normnum)
    V=pysh.expand.MakeGridGLQ(Vlm, mus, norm=normnum)
            
    return  U, V


def A22_A23(X,Y,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum):
    """Calculates equations A22 and A23 in Schwartztrauber 1996, divergence coefficients
        
        :param M: highest wavenumber
        :type M: int
        :param mus: Gaussian latitudes
        :type mus: array of float64
        
        :return: bcmn
        :rtype: array of float32
        :return: bsmn
        :rtype: array of float32
        """
    lmax = M+1
    
    Ycmnminus1 = np.zeros((lmax,lmax))
    Ycmnplus1 = np.zeros((lmax,lmax))
    Ysmnminus1 = np.zeros((lmax,lmax))
    Ysmnplus1 = np.zeros((lmax,lmax))
    brmn = np.zeros((lmax,lmax))
    bimn = np.zeros((lmax,lmax))
    
    # #Scale by 1/1-mu^2
    # X=np.multiply(X,musMAT)
    # Y=np.multiply(Y,musMAT)
    
    Xlm = pysh.expand.SHExpandGLQ(X, w, mus, norm=normnum,csphase=1,lmax_calc=lmax-1) #lmax-1?
    #print(np.shape(Xlm))
    Xcmn = np.transpose(Xlm[0,:,:])
    Xsmn = np.transpose(Xlm[1,:,:])
    
    Ylm = pysh.expand.SHExpandGLQ(Y, w, mus, norm=normnum,csphase=1,lmax_calc=lmax)  

    
    Ycmnminus1[:,1:] = np.transpose(Ylm[0,:-2,:-1])
    Ycmnplus1[:,:] = np.transpose(Ylm[0,1:,:-1])
    Ysmnminus1[:,1:] = np.transpose(Ylm[1,:-2,:-1])
    Ysmnplus1[:,:] = np.transpose(Ylm[1,1:,:-1])
    
    brmn = np.multiply(mnMAT4,Ycmnminus1) - np.multiply(mnMAT5,Ycmnplus1) - np.multiply(mnMAT1,Xsmn)
    bimn = np.multiply(mnMAT4,Ysmnminus1) - np.multiply(mnMAT5,Ysmnplus1) + np.multiply(mnMAT1,Xcmn)
                
    return brmn, bimn

def A24_A25(X,Y,M,mnMAT1,mnMAT4,mnMAT5,musMAT,w,mus,normnum):
    """Calculates equations A24 and A25 in Schwartztrauber 1996, vorticity coefficients
        
        :param M: highest wavenumber
        :type M: int
        :param mus: Gaussian latitudes
        :type mus: array of float64
        
        :return: bcmn
        :rtype: array of float32
        :return: bsmn
        :rtype: array of float32
        """
    lmax = M+1
    
    Xcmnminus1 = np.zeros((lmax,lmax))
    Xcmnplus1 = np.zeros((lmax,lmax))
    Xsmnminus1 = np.zeros((lmax,lmax))
    Xsmnplus1 = np.zeros((lmax,lmax))
    crmn = np.zeros((lmax,lmax))
    cimn = np.zeros((lmax,lmax))
    
    # #Scale by 1/1-mu^2
    # X=np.multiply(X,musMAT)
    # Y=np.multiply(Y,musMAT)
    
    Xlm = pysh.expand.SHExpandGLQ(X, w, mus, norm=normnum,csphase=1,lmax_calc=lmax)
    Xcmnminus1[:,1:] = np.transpose(Xlm[0,:-2,:-1])
    Xcmnplus1[:,:] = np.transpose(Xlm[0,1:,:-1])
    Xsmnminus1[:,1:] = np.transpose(Xlm[1,:-2,:-1])
    Xsmnplus1[:,:] = np.transpose(Xlm[1,1:,:-1])
    
    Ylm = pysh.expand.SHExpandGLQ(Y, w, mus, norm=normnum,csphase=1,lmax_calc=lmax-1)
    Ycmn = np.transpose(Ylm[0,:,:])
    Ysmn = np.transpose(Ylm[1,:,:])
    #print(np.shape(Xcmnminus1))
    
    crmn = np.multiply(mnMAT4,Xcmnminus1) - np.multiply(mnMAT5,Xcmnplus1) + np.multiply(mnMAT1,Ysmn)
    cimn = np.multiply(mnMAT4,Xsmnminus1) - np.multiply(mnMAT5,Xsmnplus1) - np.multiply(mnMAT1,Ycmn)
                
    return crmn, cimn

def f_latlon(mus,lambdas,I,J,omega,a1,test):
    fMAT=np.zeros((J,I))
    for i in range(I):
        if test==1: #Williamson Test 1
            fMAT[:,i]=mus
        elif test==2: #Williamson Test 2, so that the flow can be specified with the spherical coordinate poles no necessarily coincident with the roation axis
            for j in range(J):
                fMAT[j,i]=-np.cos(lambdas[i])*np.sqrt(1-mus[j]**2)*np.sin(a1)+(mus[j])*np.cos(a1)
    fMAT=fMAT*2*omega
    return fMAT

def step7p5(Phi1,U1,V1,w,mus,J,M,musMAT,nMAT2,normnum):
    
    PhiE=Phi1+np.multiply(0.5*(np.multiply(U1,U1)+np.multiply(V1,V1)),musMAT)
    lmax=M+1
    PhiElm = pysh.expand.SHExpandGLQ(PhiE, w, mus, norm=normnum,csphase=1,lmax_calc=lmax-1)
    
    PhiElm[0,:,:]=np.multiply(np.transpose(-nMAT2),PhiElm[0,:,:])
    PhiElm[1,:,:]=np.multiply(np.transpose(-nMAT2),PhiElm[1,:,:])
    
    PhiElmPad=np.zeros((2,J,J))
    PhiElmPad[0,:M+1,:M+1]=PhiElm[0,:,:]
    PhiElmPad[1,:M+1,:M+1]=PhiElm[1,:,:]
    
    deltaRHS2=pysh.expand.MakeGridGLQ(PhiElmPad, mus, norm=normnum)

    return deltaRHS2
    
