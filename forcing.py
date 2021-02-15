# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:28:45 2020

@author: ek672
"""

import numpy as np 


def Phieqfun(Phibar,DPhieq,lambdas,mus,I,J,g):
    PhieqMat=Phibar*np.ones((J,I)) #initialize to flat nightside temperature
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (pi,0)
            if np.pi/2<lambdas[i]<3*np.pi/2:
                PhieqMat[j,i]=PhieqMat[j,i]+DPhieq*np.cos(lambdas[i]-np.pi)*np.sqrt((1-mus[j]**2))        
    return PhieqMat


def Qfun(Phieq,Phi,taurad):
    #note Q is different from Perez-Becker and Showman, our Q is PBS-Q*g
    Q=(1/taurad)*(Phieq-Phi)
    
    return Q

def Rfun(U,V,Q,Phi, taudrag):
    Ru=np.divide(np.multiply(-U,Q),Phi)
    Rv=np.divide(np.multiply(-V,Q),Phi)
    
    #reset to 0 if losing mass
    Ru[Q<0]=0
    Rv[Q<0]=0
    
     #if taudrag is infinity, only have the R component from Perez-Becker and Showman    
    if taudrag!=-1:
        
        F=Ru-U/taudrag
        G=Rv-V/taudrag
    
    return F, G

