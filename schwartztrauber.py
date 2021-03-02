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
#import pyshtools as pysh

def mnMATgen(M,N):
    
    nMAT1=np.zeros((M+1,N+1))
    mnMAT1=np.zeros((M+1,N+1))
    mnMAT2=np.zeros((M+1,N+1))
    mnMAT3=np.zeros((M+1,N+1))
    mnMAT4=np.zeros((M+1,N+1))
    mnMAT5=np.zeros((M+1,N+1))
    
    nrange=np.arange(N+1)
    nplus1=nrange+1
    nnplus1sqrt=np.sqrt(np.multiply(nrange, nplus1))

    for m in range(M+1):
        nMAT1[m,:]=nnplus1sqrt
        
    for m in range(1,M+1):
        mnMAT1[m,1:]=m/nnplus1sqrt[1:]
        
    for m in range(M+1):
        for n in range(2,N+1):
            mnMAT2[m,n]=(n-1)*(n-m)/((2*(n-1)+1)*np.sqrt((n-1)*n))
    mnMAT2=np.triu(mnMAT2)
    
    for m in range(M+1):
        for n in range(N+1):
            mnMAT3=(n+2)*(n+1+m)/(np.sqrt((n+1)*(n+2))*(2*(n+1)+1))
            
    for m in range(M+1):
        for n in range(N+1):
            mnMAT4=(n+1)(n+m)/(np.sqrt(n*(n+1))*(2*n+1))
    
    for m in range(M+1):
        for n in range(N+1):
            mnMAT5=(n*(n-m+1))/(np.sqrt(n*(n+1))*(2*n+1))
    mnMAT5=np.triu(mnMAT2)        
    
    
    return nMAT1, mnMAT1, mnMAT2, mnMAT3

def A14(zetamn,nMAT1,mus,J):
    
    zetamnscaled=np.multiply(nMAT1,zetamn)
    zeta=pysh.expand.MakeGridGLQ(zetamnscaled, mus, lmax=J)   
    
    return zeta


def A15(deltamn,nMAT1,mus,J):
    
    deltamnscaled=np.multiply(nMAT1,-deltamn) #not sure why the minus sign
    delta=pysh.expand.MakeGridGLQ(deltamnscaled, mus, lmax=J)
    
    return delta

def A20_A21(delta,zeta,M,mnMAT1,mnMAT2,mnMAT3,w,mus,J):
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
    #temp_size = int((lmax+1)*(lmax+2)/2)
    
    # Ycmnminus1 = np.zeros((lmax+1,lmax+1))
    # Ycmnplus1 = np.zeros((lmax+1,lmax+1))
    # Ysmnminus1 = np.zeros((lmax+1,lmax+1))
    # Ysmnplus1 = np.zeros((lmax+1,lmax+1))
    # brmn = np.zeros((lmax+1,lmax+1))
    # bimn = np.zeros((lmax+1,lmax+1))
    
    deltalm = pysh.expand.SHExpandGLQ(delta, w, mus, [3,1,lmax])
    deltacmn = np.transpose(deltalm[0,:,:]) #cosine coeff
    deltasmn = np.transpose(deltalm[1,:,:]) #sine coeff
    
    zetalm = pysh.expand.SHExpandGLQ(zeta, w, mus, [3,1,lmax])
    zetacmn = np.transpose(zetalm[0,:,:])
    zetasmn = np.transpose(zetalm[1,:,:])
    
    
    
    deltacmnminus1[1:,:] = np.transpose(deltalm[0,:-2,:-1])
    deltacmnplus1[:-1,:] = np.transpose(deltalm[0,1:,:-1])
    deltasmnminus1[1:,:] = np.transpose(deltalm[1,:-2,:-1])
    deltasmnplus1[:-1,:] = np.transpose(deltalm[1,1:,:-1])
    
    zetacmnminus1[1:,:] = np.transpose(zetalm[0,:-2,:-1])
    zetacmnplus1[:-1,:] = np.transpose(zetalm[0,1:,:-1])
    zetasmnminus1[1:,:] = np.transpose(zetalm[1,:-2,:-1])
    zetasmnplus1[:-1,:] = np.transpose(zetalm[1,1:,:-1])
    
    
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
    
    ## try lmax option in makegrid, if it fails, let's pad it as we did before to satisfy the 2/3 rule
    U=pysh.expand.MakeGridGLQ(Umn, mus, lmax=J)
    V=pysh.expand.MakeGridGLQ(Vmn, mus, lmax=J)
            
    return  U, V


def A22_A23(X,Y,M,mnMAT1,mnMAT2,mnMAT3,w,mus):
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
    
    Xlm = pysh.expand.SHExpandGLQ(X, w, mus, [3,1,lmax-1])
    Xcmn = np.transpose(Xlm[0,:,:])
    Xsmn = np.transpose(Xlm[1,:,:])
    
    Ylm = pysh.expand.SHExpandGLQ(Y, w, mus, [3,1,lmax])
    Ycmnminus1[1:,:] = np.transpose(Ylm[0,:-2,:-1])
    Ycmnplus1[:-1,:] = np.transpose(Ylm[0,1:,:-1])
    Ysmnminus1[1:,:] = np.transpose(Ylm[1,:-2,:-1])
    Ysmnplus1[:-1,:] = np.transpose(Ylm[1,1:,:-1])
    
    brmn = np.multiply(mnMAT4,Ycmnminus1) - np.multiply(mnMAT5,Ycmnplus1) - np.multiply(mnMAT1,Xsmn)
    bimn = np.multiply(mnMAT4,Ysmnminus1) - np.multiply(mnMAT5,Ysmnplus1) + np.multiply(mnMAT1,Xcmn)
                
    return brmn, bimn

def A24_A25(X,Y,M,mnMAT1,mnMAT2,mnMAT3,w,mus):
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
    
    Xlm = pysh.expand.SHExpandGLQ(X, w, mus, [3,1,lmax])
    Xcmnminus1[1:,:] = np.transpose(Xlm[0,:-2,:-1])
    Xcmnplus1[:-1,:] = np.transpose(Xlm[0,1:,:-1])
    Xsmnminus1[1:,:] = np.transpose(Xlm[1,:-2,:-1])
    Xsmnplus1[:-1,:] = np.transpose(Xlm[1,1:,:-1])
    
    Ylm = pysh.expand.SHExpandGLQ(Y, w, mus, [3,1,lmax-1])
    Ycmn = np.transpose(Ylm[0,:,:])
    Ysmn = np.transpose(Ylm[1,:,:])
    
    crmn = np.multiply(mnMAT4,Xcmnminus1) - np.multiply(mnMAT5,Xcmnplus1) + np.multiply(mnMAT1,Ysmn)
    cimn = np.multiply(mnMAT4,Xsmnminus1) - np.multiply(mnMAT5,Xsmnplus1) - np.multiply(mnMAT1,Ycmn)
                
    return crmn, cimn
