# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:37:55 2021

@author: ek672
"""
import numpy as np 

y=np.arange(5)

print(y[4])
print(y[:4])

x=np.zeros(5)

print(len(y[:3]))
x[:3]=y[:3]
print(x)


M=63

G=np.ones((M+1,M+1))

for m in range(0,M+1):
    for n in range(1,M+1):
        G[m,n]=np.arctan(M/2-m)/np.pi+0.5
        
import matplotlib.pyplot as plt
import metview as mv

plt.plot (np.arange(M+1),G[:,2])
