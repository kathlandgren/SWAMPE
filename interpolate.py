# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 14:29:07 2022

@author: ek672
"""

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt




points=[[1,0.1],[1,1],[1,10],[5,0.1],[5,1],[5,10],[10,0.1],[10,1],[10,10]]


values_high=[0.62,	1.60,	6.03,
0.36,	2.38,	3.75,
0.23,	1.45,	5.87]

values_med=[0.42,	0.73,	4.27,
0.31,	2.11,	3.66,
0.18,	1.24,	7.83]

values_low=[0.11,	0.20,	2.19,
0.16,	0.36,	4.66,
0.13,	0.88,	4.59]

grid_x, grid_y = np.mgrid[1:10:1000j, 0.1:10:1000j]

met='linear'

grid_high= griddata(points, values_high, (grid_x, grid_y), method=met)

grid_med = griddata(points, values_med, (grid_x, grid_y), method=met)

grid_low = griddata(points, values_low, (grid_x, grid_y), method=met)

# plt.contour(grid_x, grid_y, grid_high, levels=[1],colors=('k',),linestyles=('-',),linewidths=(2,))

# plt.contour(grid_x, grid_y, grid_med, levels=[1],colors=('k',),linestyles=('--',),linewidths=(2,))

# plt.contour(grid_x, grid_y, grid_low, levels=[1],colors=('k',),linestyles=('-.',),linewidths=(2,))

# plt.xlabel(r"Rotation period $P_{rot}$")

# plt.ylabel(r"Radiative timescale $\tau_{{rad}}$")

# plt.scatter([1,1,1,5,5,5,10,10,10],3*[0.1,1,10],marker='x')

# plt.show()
fig=plt.figure(figsize=((6,6)),dpi=300)
CS0=plt.contourf((grid_x), (grid_y))
CS1=plt.contour((grid_x), (grid_y), (grid_high), levels=[1], colors=('k',),linestyles=('-',),linewidths=(2,))
CS0p5=plt.contour((grid_x), (grid_y), (grid_med), levels=[1],colors=('k',),linestyles=('--',),linewidths=(2,))
CS0p1=plt.contour((grid_x), (grid_y), (grid_low), levels=[1],colors=('k',),linestyles=('-.',),linewidths=(2,))

manual_locations_high = [(2,1)]
manual_locations_med = [(3,1.5)]
manual_locations_low = [(3,1)]
plt.clabel(CS1, inline=1, fontsize=8,inline_spacing=-30,use_clabeltext=1,fmt=r'$\Delta \Phi_{\rm eq}/\overline{\Phi}=1$',manual=manual_locations_high)
plt.clabel(CS0p5, inline=1, fontsize=8,inline_spacing=-10,fmt=r'$\Delta \Phi_{\rm eq}/\overline{\Phi}=0.5$',manual=manual_locations_med)
plt.clabel(CS0p1, inline=1, fontsize=8,inline_spacing=0,fmt=r'$\Delta \Phi_{\rm eq}/\overline{\Phi}=0.1$',manual=manual_locations_low)
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.05,12)
plt.xlim(0.95,12)
plt.xlabel(r"Rotation period $P_{rot}$ [days]")

plt.ylabel(r"Radiative timescale $\tau_{{rad}}$ [days]")
plt.title(r"Contours where $\tau_{\rm rad}=\tau_{\rm adv}$")



plt.scatter([1,1,1,5,5,5,10,10,10],3*[0.1,1,10],marker='x')

plt.show()
#plt.savefig('tau_adv_contours.pdf')
#plt.scatter(np.log10([1,1,1,5,5,5,10,10,10]),np.log10([0.1,1,10, 0.1,1,10,0.1,1,10]))

