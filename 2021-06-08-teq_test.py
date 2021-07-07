# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 14:46:44 2021

@author: ek672
"""

## imports
import initial_conditions as ic
import forcing
import params as p
import testing_plots
import numpy as np
from astropy import constants as const
import continuation as cont

## define necessary dimensiony things

M=42

N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

#physical params
g=p.g
Phibar=p.Phibar
DPhieq=Phibar

#DGDPhieq=Phibar#2.5*Phibar

R=p.R
k1=0.0002
k2=0.0004
pressure=p.pressure
Cp=p.Cp
sigmaSB=p.sigmaSB

taurad=p.taurad

# distance=0.05*const.au.to_value()
#L_star= (10**0.25)*(const.L_sun.to_value()) # HD 209458
#DGPhibar=1450*R #Teq from Sing et al from exoplanet archive

distance=0.03*const.au.to_value()#0.02*const.au.to_value() # WASP 18
L_star=  1.29*(10**26)#(10**0.39)*(const.L_sun.to_value()) 
DGPhibar=1000*3700#Phibar# 2400*R #Teq from Salz et al 2015 from exoplanet archive


DGDPhiEq=((L_star)/(4*np.pi*(distance**2)))/sigmaSB #2000000# 924637 

def insolation(Lstar,distance,sigmaSB):
    insolation=((L_star)/(4*np.pi*(distance**2)))/sigmaSB
    
    return insolation

## load Phidata
tindex=68
Phidata = cont.load_and_restore('data\Showman2015\Phidata-k1-0.0002-k2-0.0004',I)
Phiic0=Phidata[tindex,:,:]


## PBS Teq
PhiEq=forcing.Phieqfun(Phibar,DPhieq,lambdas,mus,I,J,g)
PBSTeq=testing_plots.geopot_to_temp(PhiEq, R)

PBSPhi_forcing=forcing.Qfun(PhiEq,Phiic0,Phibar,taurad)

## Langton Teq

DGTeq_exp=forcing.DoubleGrayTEqfun(DGPhibar,DGDPhiEq,lambdas,mus,I,J,k1,k2,pressure,g,R,Cp,sigmaSB)
DGTeq=DGTeq_exp**0.25

DGPhi_forcing=forcing.DoubleGrayPhiForcing(DGTeq_exp,Phiic0,Phibar,k2,sigmaSB,Cp,R)

##Plot
test=p.test
minlevel=1000
maxlevel=2600
t=0
dt=1
a1=p.a1

# testing_plots.temp_plot(PBSTeq,lambdas,mus,t,dt,10,a1,minlevel,maxlevel)
# testing_plots.temp_plot(DGTeq,lambdas,mus,t,dt,test,a1,700,1600)
# testing_plots.temp_plot(np.abs(PBSTeq-DGTeq),lambdas,mus,t,dt,test,a1,0,400)
# testing_plots.physical_plot(PBSTeq, mus, lambdas)
# testing_plots.physical_plot(DGTeq, mus, lambdas)

#testing_plots.physical_plot(DGTeq, mus, lambdas)
testing_plots.physical_plot(PBSPhi_forcing, mus, lambdas)
testing_plots.physical_plot(DGPhi_forcing, mus, lambdas)
testing_plots.physical_plot(Phiic0, mus, lambdas)

testing_plots.temp_plot(PBSPhi_forcing,lambdas,mus,t,dt,10,a1,0,50)
   

