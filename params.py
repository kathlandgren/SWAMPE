
"""
Created on Wed Jun  3 18:41:45 2020

@author: ek672

This is the file containing the spectral, physical, and code parameters.
"""

import numpy as np

# #Spectral parameters
M=63 #the largest Fourier wave number
# N=M #highest degree of the Legendre functions for m=0
# I=192 #length of array/ number of samples
# J=96#int(np.ceil(I/2))


#time-stepping parameters
tmax=200 #number of time steps
# dt=900 #time step length, in seconds


#make these into a file that gets read later
test=10
##specifies the testing regime: 
# 1 -- test 1 from Williamson (advection of cosine bell)
# 2 -- test 2 from Williamson (global steady state nonlinear zonal geostrophic flow)
# 3 -- test 3 from Williamson (global steady state nonlinear zonal geostrophic flow with compact support)
# 6 -- test 6 from Williamson (Rossby-Haurwitz wave)
# 10 -- Hot Jupiter (PBS) 
a1=np.pi/2 #alpha from Test 1 and 2

if test==1: # Williamson Test 1
    forcflag=0
    omega=7.2921159*10**(-5) #rotation rate of the planet, radians per second
    a=6.37122*10**(6)  #radius of the planet, meters
    Phibar=1*(10**3) #Geopotential height
    g=9.8 #gravity of the planet in m/s^2
    
    minlevel=3 #the log values for the colorbar plotting.
    maxlevel=3.3
elif test==2: # Williamson Test 2
    forcflag=0
    omega=7.2921159*10**(-5) #rotation rate of the planet, radians per second
    a=6.37122*10**(6)  #radius of the planet, meters
    Phibar=3*(10**3) #Geopotential height m
    g=9.8 #gravity of the planet in m/s^2
    
    minlevel=3 #the log values for the colorbar plotting.
    maxlevel=4
elif test==10: # PBS Hot Jupiter
    #Physical parameters
    forcflag=1
    omega=3.2*(10**(-5))#7.2921159*10**(-5) #rotation rate of the planet, radians per second
    a=8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
    Phibar=4*(10**6) #1*(10**3) #Geopotential height #maybe 10*6 instead>
    g=9.8
    
    minlevel=4 #np.log10(Phibar) should be good #the log values for the colorbar plotting.
    maxlevel=7
#Hyperviscosity parameters
diffflag=1

#Modal Splitting Fiter 
modalflag=1
alpha=0.01 #filter coefficient to prevent aliasing

#forcing parameters
taurad=3600*24*100#in Earth days
taudrag=3600*24*100#-1 #3600*24*1 #if set to -1, means infinity
DPhieq=Phibar


expflag=1
#1 means explicit,
#anything else means semi-implicit scheme

#lat-lon grid points
# lambdas=np.linspace(-np.pi, np.pi-1/I, num=I) #longitudes
# [mus,w]=sp.roots_legendre(J) #Gaussian latitudes and weights
#\mu ranges from -1 to 1,\lambda ranges from 0 to 2\pi

