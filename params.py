
"""
Created on Wed Jun  3 18:41:45 2020

@author: ek672

This is the file containing the spectral, physical, and code parameters.
"""

import numpy as np

# #Spectral parameters
M=63 #the largest Fourier wave number

#time-stepping parameters

tmax=500#5000#864 #number of time steps

# dt=900 #time step length, in seconds


# #make these into a file that gets read later
# test=1
# ##specifies the testing regime: 
# #1 -- test 1 (advection)
# # 2 -- Hot Jupiter (PBS) 
# a1=np.pi/2-0.05 #alpha from test 1

# if test==1: 
#     omega=7.2921159*10**(-5) #rotation rate of the planet, radians per second
#     a=6.37122*10**(6)  #radius of the planet, meters
#     Phibar=1*(10**3) #Geopotential height
#     g=9.8 #gravity of the planet in m/s^2
# elif test==2:
#     #Physical parameters
#     omega=3.2*(10**(-5))#7.2921159*10**(-5) #rotation rate of the planet, radians per second
#     a=8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
#     Phibar=4*(10**6) #1*(10**3) #Geopotential height #maybe 10*6 instead>
#     g=9.8


#make these into a file that gets read later
test=11
##specifies the testing regime: 
# 1 -- test 1 from Williamson (advection of cosine bell)
# 2 -- test 2 from Williamson (global steady state nonlinear zonal geostrophic flow)
# 3 -- test 3 from Williamson (global steady state nonlinear zonal geostrophic flow with compact support)
# 6 -- test 6 from Williamson (Rossby-Haurwitz wave)
# 10 -- Hot Jupiter (PBS) 
a1=0.0#np.pi/2 #alpha from Test 1 and 2

if test==1: # Williamson Test 1
    forcflag=0
    expflag=1 #1 means explicit, anything else means semi-implicit scheme
    omega=7.2921159*10**(-5) #rotation rate of the planet, radians per second
    a=6.37122*10**(6)  #radius of the planet, meters
    Phibar=1*(10**3) #Geopotential height
    g=9.8 #gravity of the planet in m/s^2
    
    minlevel=3 #the log values for the colorbar plotting.
    maxlevel=3.3
elif test==2: # Williamson Test 2
    forcflag=0
    expflag=1 #1 means explicit, anything else means semi-implicit scheme
    omega=7.2921159*10**(-5) #rotation rate of the planet, radians per second
    a=6.37122*10**(6)  #radius of the planet, meters
    Phibar=3*(10**3) #Geopotential height m
    g=9.8 #gravity of the planet in m/s^2
    
    minlevel=3 #the log values for the colorbar plotting.
    maxlevel=5
elif test==10: # PBS Hot Jupiter
    #Physical parameters
    forcflag=1
    expflag=0 #1 means explicit, anything else means semi-implicit scheme
    omega=3.2*10**(-5) #1.46*10**(-5) #rotation rate of the planet, radians per second
    a=8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
    Phibar=4*(10**6) #1*(10**3) #Geopotential height #maybe 10*6 instead>
    g=9.8
    
    minlevel=6.55 #np.log10(Phibar) should be good #the log values for the colorbar plotting.
    maxlevel=6.8
    
elif test==11: #Langton hot Jupiter
    #Physical parameters
    forcflag=1
    expflag=0 #1 means explicit, anything else means semi-implicit scheme
    omega=3.2*10**(-5) #1.46*10**(-5) #rotation rate of the planet, radians per second
    a=8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
    Phibar=4*(10**6) #1*(10**3) #Geopotential height #maybe 10*6 instead>
    DPhieq=Phibar
    g=9.8
    k1=2*10**(-4)
    k2=4*10**(-4)
    pressure=100*250*g/10 #(in Pa)
    R=3000
    Cp=13000 
    sigmaSB=5.7*10**(-8)
    
    minlevel=6.55 #np.log10(Phibar) should be good #the log values for the colorbar plotting.
    maxlevel=6.8

#Continuation flag to load
contflag=0  #continuation flag to save
saveflag=0
#Continuation save frequency: every savefreq time steps
savefreq=100

#Hyperviscosity parameters
diffflag=1

#Modal Splitting Fiter 
modalflag=1
alpha=0.01 #filter coefficient to prevent aliasing

#Plotting flag
plotflag=1
#plotting frequency, every plotfreq frames
plotfreq=5

#forcing parameters

taurad=3600*24*10 #in Earth days
taudrag=-1#3600*24*10#if set to -1, means infinity


DPhieq=Phibar

