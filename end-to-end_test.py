# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 15:33:22 2022

@author: ek672

This script runs end-to-end test
"""
 
import main_function_restructure as main
import numpy as np
import continuation as cont

### TEST 1 ###

#spectral resolution
M=42

dt=1200
# Length of the run in time steps
tmax = 1000#20

#mean geopotential height. In hot Jupiter case, Phibar is the flat nightside thickness
Phibar=1000

#rotation rate of the planet, radians per second
omega=7.2921159*10**(-5)
#planetary radius, meters
a=6.37122*10**(6)
#angle for test cases 1 and 2, radians
a1=0.0


# main.main(M,dt,tmax,Phibar, omega, a, a1=0.0, test=1,plotflag=1, forcflag=0, plotfreq=10, saveflag=0, minlevel=900,maxlevel=1800,diffflag=0,modalflag=0)

# #save this somehow? -- low save frequency! 

# #load reference data?

# #assert the data are close? 

# main.main(M,dt,tmax,Phibar, omega, a, a1=np.pi/3, test=1,plotflag=1, forcflag=0, plotfreq=10, saveflag=0, minlevel=900,maxlevel=1800,diffflag=0,modalflag=0)



# ### TEST 2 ###
# M=63
# Phibar=3000

# main.main(M,dt,tmax,Phibar, omega, a, a1=0.05, test=2,plotflag=1, plotfreq=10,forcflag=0, saveflag=0, minlevel=-10000,maxlevel=10**4,diffflag=1,modalflag=1)

# main.main(M,dt,tmax,Phibar, omega, a, a1=np.pi/4, test=2,plotflag=1, plotfreq=10,forcflag=0, saveflag=0, minlevel=-10000,maxlevel=10**4,diffflag=1,modalflag=1)


### TEST TERRESTRIAL PLANET ###
main.main(M,dt,tmax,Phibar, omega, a, a1=np.pi/4, test=9,plotflag=1, plotfreq=10,forcflag=1, DPhieq=1500, saveflag=0, minlevel=-10000,maxlevel=10**4,diffflag=0,modalflag=1)

### TEST HOT JUPITER ###

M=42
dt=180
#paste hot jupiter test values in here
#main.main(M,dt1,tmax,Phibar, omega, a, test=p.test, DPhieq=p.DPhieq, plotflag=p.plotflag, forcflag=p.forcflag, plotfreq=p.plotfreq, minlevel=p.minlevel, maxlevel=p.maxlevel, saveflag=p.saveflag, savefreq=p.savefreq, k1=k1, k2=k2,taudrag=p.taudrag, taurad=p.taurad, diffflag=p.diffflag,g=p.g,alpha=p.alpha,K6=K6)