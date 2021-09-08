
"""
Created on Wed Jun  3 18:41:45 2020

@author: ek672

This is the file containing the spectral, physical, and code parameters.
"""

import numpy as np

# #Spectral parameters
M=42 #the largest Fourier wave number

#time-stepping parameters

tmax=1000#5000#864 #number of time steps

# dt=900 #time step length, in seconds

#make these into a file that gets read later
test=9
##specifies the testing regime: 
# 1 -- test 1 from Williamson (advection of cosine bell)
# 2 -- test 2 from Williamson (global steady state nonlinear zonal geostrophic flow)
# 3 -- test 3 from Williamson (global steady state nonlinear zonal geostrophic flow with compact support)
# 6 -- test 6 from Williamson (Rossby-Haurwitz wave)
# 10 -- Hot Jupiter (PBS) 
# 11 -- Double-gray scheme
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
    
elif test==9: # Kraucunas and Hartman 2006 basic state
    #Physical parameters
    forcflag=1
    expflag=0 #1 means explicit, anything else means semi-implicit scheme
    omega=7.2921159*10**(-5)#3.2*10**(-5) #1.46*10**(-5) #rotation rate of the planet, radians per second
    a=6.37122*10**(6)#8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
    Phibar=1#4*(10**6) #1*(10**3) #Geopotential height 
    g=9.8 #m/s
    DPhieq=1500#4*(10**6) #m^2/s^2
    
    minlevel=6.55 #np.log10(Phibar) should be good #the log values for the colorbar plotting.
    maxlevel=6.8
    
elif test==10: # PBS Hot Jupiter
    #Physical parameters
    forcflag=1
    expflag=0 #1 means explicit, anything else means semi-implicit scheme
    omega=3.2*10**(-5) #1.46*10**(-5) #rotation rate of the planet, radians per second
    a=8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
    Phibar=4*(10**6) #1*(10**3) #Geopotential height 
    g=9.8 #m/s
    DPhieq=0.001*Phibar#4*(10**6) #m^2/s^2
    
    minlevel=6.55 #np.log10(Phibar) should be good #the log values for the colorbar plotting.
    maxlevel=6.8
elif test==11: #Langton hot Jupiter -- for DOuble Gray Forcing
    #Physical parameters
    forcflag=1
    expflag=0 #1 means explicit, anything else means semi-implicit scheme
    omega=3.3*10**(-5)#8.25*10**(-6)#3.3*10**(-5)#3.2*10**(-5) #1.46*10**(-5) #rotation rate of the planet, radians per second
    a=8.2*(10**7)#6.37122*10**(6)  #radius of the planet, meters
    Phibar=1000*3700#4*(10**6) #Geopotential height based on T=500K (Phi=RT) #1.5*(10**6)
    DPhieq=8941526274934.578#1257402132412##Phibar#2.5*Phibar#2*Phibar#2.31*(10**6) #Estimated at 0.3 AU from the Sun #Phibar 

    g=9.8
    k1=2*10**(-4) #2 is used in Langton and Laughlin
    k2=4*10**(-4)
    pressure=100*250*g/10 #(in Pa)
    R=3700#3000
    Cp=13000 
    sigmaSB=5.7*10**(-8)
    
    minlevel=6.55#np.log10(Phibar) should be good #the log values for the colorbar plotting.
    maxlevel=6.8

#Continuation flag to load
contflag=0 
#continuation flag to save
saveflag=0
#Continuation save frequency: every savefreq time steps
savefreq=180

#Hyperviscosity parameters
diffflag=1

#Modal Splitting Fiter 
modalflag=1
alpha=0.01 #filter coefficient to prevent aliasing

#Plotting flag
plotflag=1
#plotting frequency, every plotfreq frames
plotfreq=50

#forcing parameters

taurad=int(3600*24*10) #in Earth days
taudrag=int(3600*24*10)#if set to -1, means infinity
