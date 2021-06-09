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


## define necessary dimensiony things

M=63

N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)

#physical params
g=p.g
Phibar=p.Phibar
DPhieq=Phibar

DGDPhieq=Phibar#2.5*Phibar

R=p.R
k1=0.0002
k2=0.0004
pressure=p.pressure
Cp=p.Cp
sigmaSB=p.sigmaSB


## PBS Teq
PhiEq=forcing.Phieqfun(Phibar,DPhieq,lambdas,mus,I,J,g)
PBSTeq=testing_plots.geopot_to_temp(PhiEq, R)



## Langton Teq

DGTeq_exp=forcing.DoubleGrayTEqfun(Phibar,DGDPhieq,lambdas,mus,I,J,k1,k2,pressure,g,R,Cp,sigmaSB)
DGTeq=DGTeq_exp**0.25

##Plot
test=p.test
minlevel=1000
maxlevel=2600
t=0
dt=1
a1=p.a1

testing_plots.temp_plot(PBSTeq,lambdas,mus,t,dt,10,a1,minlevel,maxlevel)
testing_plots.temp_plot(DGTeq,lambdas,mus,t,dt,test,a1,minlevel,maxlevel)
testing_plots.temp_plot(PBSTeq-DGTeq,lambdas,mus,t,dt,test,a1,0,400)
# testing_plots.physical_plot(PBSTeq, mus, lambdas)
# testing_plots.physical_plot(DGTeq, mus, lambdas)

testing_plots.physical_plot(PBSTeq-DGTeq, mus, lambdas)
   

