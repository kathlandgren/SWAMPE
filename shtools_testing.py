# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 09:38:05 2021

@author: ek672
"""
import initial_conditions as ic
import numpy as np
import params as p
import pyshtools as pysh
import testing_plots


M=63
a=p.a
omega=p.omega
a1=np.pi/2

normnum=1
samplingnum=2

lmax=M

N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)


SU0, sina, cosa, etaamp, Phiamp =ic.test1_init(a, omega, a1)

Uic=np.full((J,I),0.0) #initialize
Vic=np.full((J,I),0.0)
    

for i in range(I):
    for j in range(J):
        Uic[j,i]=SU0*(np.cos(np.arcsin(mus[j]))*cosa +mus[j]*np.cos(lambdas[i])*sina)#assign according to init.f
        Vic[j,i]=-SU0*np.sin(lambdas[i])*sina



Umn=pysh.expand.SHExpandGLQ(Uic, w, mus,norm=normnum,lmax_calc=lmax)
Vmn=pysh.expand.SHExpandGLQ(Vic, w, mus,norm=normnum,lmax_calc=lmax)

UmnPad=np.zeros((2,J,J)) #padding so that shtools give us the correct dimensions for lat lon
UmnPad[:,:M+1,:M+1]=Umn

VmnPad=np.zeros((2,J,J)) #padding so that shtools give us the correct dimensions for lat lon
VmnPad[:,:M+1,:M+1]=Vmn


Unew=pysh.expand.MakeGridGLQ(UmnPad, mus, norm=normnum)
Vnew=pysh.expand.MakeGridGLQ(VmnPad, mus, norm=normnum)

testing_plots.physical_plot(Uic-Unew,mus,lambdas)
#testing_plots.physical_plot(Vic-Vnew,mus,lambdas)
print(np.max(np.abs(Uic-Unew)))


#### Lat-lon SHtools

#define latitudes

latfull=np.linspace(np.pi/2,-np.pi/2,num=J+1,endpoint=True)
lat=latfull[:J]

lonfull=np.linspace(0.0,2*np.pi,num=I+2,endpoint=True)
lon=lonfull[:I+1]


#Uiclatlon=np.full((J,J),0.0) #initialize

Uiclatlon=np.full((J,I+1),0.0) #initialize
Viclatlon=np.full((J,I+1),0.0)
    

for i in range(I+1):
    for j in range(J):
        Uiclatlon[j,i]=SU0*(np.cos(lat[j])*cosa +np.sin(lat[j])*np.cos(lon[i])*sina)#assign according to init.f
        Viclatlon[j,i]=-SU0*np.sin(lon[i])*sina



Umnlatlon = pysh.expand.SHExpandDH(Uiclatlon, norm=normnum, sampling=samplingnum)# lmax_calc=lmax)


Ulatlonnew = pysh.expand.MakeGridDH(Umnlatlon,norm=normnum, sampling=samplingnum)

testing_plots.physical_plot_latlon(Uiclatlon-Ulatlonnew,lat,lon)

print(np.max(np.abs(Uiclatlon-Ulatlonnew)))