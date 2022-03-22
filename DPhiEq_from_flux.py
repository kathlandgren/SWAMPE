# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 19:20:00 2022

@author: ek672
"""
import params as p

flux1=468000
flux2=7.37*10**4
flux3=1.16*10**4

DPhiEq1=3700*(flux1/(4*5.7*10**(-8)))**0.25/p.Phibar
DPhiEq2=3700*(flux2/(4*5.7*10**(-8)))**0.25/p.Phibar
DPhiEq3=3700*(flux3/(4*5.7*10**(-8)))**0.25/p.Phibar



g=21.1
p=100*250#*g/10 #(in Pa)
Cp=13000

sigmaSB=5.7*10**(-8)

#equilibrium temp
TeH=1198
TeW=755

def find_tau_rad(p,Cp,g,sigmaSB,Te):
    
    taurad=((p*Cp)/(4*g*sigmaSB*(Te**3)))/(3600*24)
    
    return taurad

tauradH=find_tau_rad(p, Cp, g, sigmaSB, TeH)
tauradW=find_tau_rad(p, Cp, g, sigmaSB, TeW)