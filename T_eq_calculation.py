# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:03:41 2022

@author: ek672
"""

import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt


#implement Teq temperature

def get_Teq(A,Sp):
    Teq4=((1-A)*Sp/4)/const.sigma_sb.value
    Teq=(Teq4)**0.25
    return Teq

TeqEarth=get_Teq(0.3,1361)
#print(TeqEarth)

def get_Sp(L,a):
    Sp=L/(4*np.pi*a**2)
    return Sp
    
SpEarth=get_Sp(const.L_sun.value,const.au.value)
#print(SpEarth)

#implement third law for every body

def get_rotation_period(a,M):
    
    psquared=((4*np.pi**2)/(const.G.value*M))*a**3
    p=np.sqrt(psquared)
    
    return p

def get_distance(p,M):
    
    acubed=(p**2)*(const.G.value*M)/(4*np.pi**2)
    a=acubed**(1/3)
    
    return a

#test for Earth
p=get_rotation_period(const.au.value, const.M_sun.value)
#print(p/(3600*24))

a=get_distance(31558196.02038122, const.M_sun.value)
#print(a/const.au.value)

def get_luminosity(Rstar,Tstar):
    L=4*np.pi*Rstar**2*const.sigma_sb.value*Tstar**4

def get_Teq_and_dist(startype,prot,albedo):
    
    if startype=='Sun': 
        Mstar=const.M_sun.value
        Lstar=const.L_sun.value   
    elif startype=='K0': 
        Mstar=0.78*const.M_sun.value
        Tstar=5240
        Rstar=0.85*const.R_sun.value
        Lstar=0.40*const.L_sun.value    
    elif startype=='K5': 
        Mstar=0.69*const.M_sun.value
        Tstar=4410
        Rstar=0.74*const.R_sun.value
        Lstar=0.16*const.L_sun.value
    elif startype=='M0':
        Mstar=0.60*const.M_sun.value
        Tstar=3800
        Rstar=0.51*const.R_sun.value
        Lstar=0.072*const.L_sun.value
    elif startype=='M5': 
        Mstar=0.15*const.M_sun.value
        Tstar=3120
        Rstar=0.18*const.R_sun.value
        Lstar=0.0027*const.L_sun.value
    else:
        print('Unsupported star type.')
        
        #get distance
    a=get_distance(prot,Mstar)
    #print('distance is '+str(a/const.au.value))
    Sp=get_Sp(Lstar,a)
    Teq=get_Teq(albedo,Sp)
    
    distance=a/const.au.value
    # print('distance is '+str(a/const.au.value)+' AU')
    # print('temperature is '+str(Teq))
    
    return Teq, distance

prot=10

typelist=['K0','K5', 'M0', 'M5']
protlist=[24*3600, 5*24*3600, 10*24*3600]
Teqlist=np.zeros((len(typelist),len(protlist)))
distlist=np.zeros((len(typelist),len(protlist)))

for i in range(len(typelist)):
    for j in range(len(protlist)):
        Teq,distance=get_Teq_and_dist(typelist[i], protlist[j], 0)
        Teqlist[i,j]=Teq
        distlist[i,j]=distance
        
print(Teqlist[0,:])

protlisth=[24, 120, 240]

plt.scatter(protlisth,Teqlist[0],c='purple')
plt.scatter(protlisth,Teqlist[1],c='blue')
plt.scatter(protlisth,Teqlist[2],c='orange')
plt.scatter(protlisth,Teqlist[3],c='red')
plt.xlabel(r'Planetary rotation/orbital period $P_{rot}}$, hours')
plt.ylabel(r"Equilibrium temperature $T_{eq}$, K")
plt.grid(visible=True, which='major', axis='both')
plt.fill_between(protlisth,400, 1200, facecolor='gray', interpolate=True, alpha=0.5)
plt.text(150, 1000, 'sampled region', style='italic',
        bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 10})
#ax.set_yticks([0.2, 0.6, 0.8], minor=False)
plt.xticks([24, 120, 240])
plt.yticks([400, 800, 1200, 1600, 2000])
plt.legend(typelist)
plt.savefig('t_eq.pdf', dpi = 300,bbox_inches='tight')
plt.show()