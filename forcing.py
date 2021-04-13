# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:28:45 2020

@author: ek672
"""

import numpy as np 

# def heqfun(Phibar,Dheq,lambdas,mus,I,J,g):
#     heqMat=Phibar*np.ones((J,I))/g #H
#     #heqMat=np.zeros((J,I))
    
#     for i in range(I):
#         for j in range(J):
#             #assume substellar point is (0,0)
#             if -np.pi/2<lambdas[i]<np.pi/2:
#                 heqMat[j,i]=heqMat[j,i]+Dheq*np.cos(lambdas[i])*(1-mus[j]**2)
            
#     return heqMat



def Phieqfun(Phibar,DPhieq,lambdas,mus,I,J,g):

    PhieqMat=Phibar*np.ones((J,I)) #initialize to flat nightside temperature
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (0,0)
            if  -np.pi/2<lambdas[i]<np.pi/2:
                PhieqMat[j,i]=PhieqMat[j,i]+DPhieq*np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))        
    return PhieqMat



# def Qfun(heq,Phi,Phibar,taurad,g):
#     Q=(1/taurad)*(heq-(Phi+Phibar)/g)
#     #Q=(1/taurad)*(heq-(Phi)/g)
#     return Q
def DoubleGrayPhiEqfun(Phibar,DPhieq,lambdas,mus,I,J,k1,k2,p,g,R,Cp,sigma):
    
    TeqMat=((Phibar/R)**4)*np.ones((J,I))

    x=np.exp(-k1*p/g)
    print(x)
    
    for i in range(I):
        for j in range(J):
            #assume substellar point is (0,0)
            if  -np.pi/2<lambdas[i]<np.pi/2:
                
                TeqMat[j,i]=((k1/k2)*((DPhieq+Phibar)/R)**4*x**(1/(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))))+(Phibar/R)**4)
                #PhieqMat[j,i]=PhieqMat[j,i]+k1*DPhieq*(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2)))*x**(1/(np.cos(lambdas[i])*np.sqrt((1-mus[j]**2))))        
    return TeqMat

def DoubleGrayPhiForcing(TeqMat,Phidata,Phibar,k2,sigma,Cp,R):
    outer_coeff=sigma*k2*R/Cp
    Q=outer_coeff*(TeqMat-((Phidata+Phibar)/R)**4)
    
    return Q
    
    

def Qfun(Phieq,Phi,Phibar,taurad):
    #note Q is different from Perez-Becker and Showman, our Q is PBS-Q*g
    Q=(1/taurad)*(Phieq-(Phi+Phibar))
    
    return Q


def Qfun_with_rampup(Phieq,Phi,Phibar,taurad,t,dt):
    #note Q is different from Perez-Becker and Showman, our Q is PBS-Q*g
    
    #slowly ramp up
    if t*dt<15*3600:
        factor=t*dt/(15*3600)
    else:
        factor=1
    

        
    Q=factor*(1/taurad)*(Phieq-(Phi+Phibar))
    print(np.max(Q))
    
    return Q


def Rfun(U,V,Q,Phi,Phibar, taudrag):
    Ru=np.divide(np.multiply(-U,Q),Phi+Phibar)
    Rv=np.divide(np.multiply(-V,Q),Phi+Phibar)
    
    #reset to 0 if losing mass
    Ru[Q<0]=0
    Rv[Q<0]=0
    
     #if taudrag is infinity, only have the R component from Perez-Becker and Showman    
    if taudrag!=-1:
        
        F=Ru-(U/taudrag)
        G=Rv-(V/taudrag)
        
    else:
        F=Ru
        G=Rv
        
    
    return F, G
