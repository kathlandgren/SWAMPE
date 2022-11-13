# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:44:57 2022

@author: ezets
"""

import numpy
import SWAMPE
import time
start_time = time.time()


# planetary radius a in meters

a=8.2*10**(7) 

# planetary rotation rate omega in radians/s

omega=3.2*10**(-5)

# the reference geopotential height, m^2/s^2, typically based on scale height

Phibar=4*(10**6)

#timestep length in seconds
dt=30
#number of timesteps for the simulation
tmax=180000

M=42

DPhieq=Phibar

taurad=int(3600*24) 
taudrag=-1

plotflag=False

SWAMPE.run_model(M,dt,tmax,Phibar, omega, a,  taurad=taurad, taudrag=taudrag, DPhieq=DPhieq, plotflag=plotflag, saveflag=True,savefreq=1200,verbose=True,timeunits='hours',K6=1.24*10**35) 

print("--- %s seconds ---" % (time.time() - start_time))