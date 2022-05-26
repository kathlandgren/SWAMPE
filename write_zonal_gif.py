# -*- coding: utf-8 -*-
"""
Created on Mon May 23 18:01:06 2022

@author: ek672
"""


import matplotlib.pyplot as plt
import imageio
import numpy as np


import initial_conditions as ic
import pickle

M=42
#get other dimensional parameters using the spectral dimension
N,I,J,dt,K4,lambdas,mus,w=ic.spectral_params(M)
#Pmn, Hmn = rfl.PmnHmn(J, M, N, mus)
        
forcingarray=['high','medium','low']
phibararray=['1','2','3','4']
tauarray=['0p1','1','10']
periodarray=['1','5','10']

tmax=4780
tmin=2780

assert tmax>tmin, "tmax must be greater than tmin"

trange=np.arange(tmin,tmax)
spacing=10
tspaced=trange[::spacing]



def read_pickle(path,filename):
    infile = open(path+filename,'rb')
    var = pickle.load(infile)
    infile.close()
    return var


def write_pickle(filename,data):
    outfile = open(filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()



num_snapshot=len(tspaced)


def write_zonal_gif(mus,Phidata,tmax,frms,string):
    """ Writes a .gif file using the plot_for_offset function
    :param z_max: max plot height
    :type z_max: float
    :param lambdas: longitudes
    :type lambdas: list
    :param mus: Gaussian latitudes
    :type mus: list
    :param Xidata: state variable data for all time
    :type Xidata: list
    :param tmax: maximal time
    :type tmax: int
    :param frms: frames per second
    :type frms: int
    :param string: gif name
    :type string: string
    """
    #kwargs_write = {'fps':1.0, 'quantizer':'nq'}
    imageio.mimsave('./'+string, [zonal_wind_plot(Phidata[i,:,:],mus,i) for i in range(tmax)], fps=frms)



    
def zonal_wind_plot(plotdata,mus,t):
    
    """

    :param plotdata: first variable, JxM+1 
    :type statevar1: float

    
    """
    #t = np.linspace(0, dt*tmax/3600, tmax, endpoint=True)
    fig,ax= plt.subplots(figsize=(5,5))
    Ubar=np.mean(plotdata,axis=1)

    Y = np.arcsin(mus)*180/np.pi #generate latitude
    
    plt.plot(Ubar, Y)
    
    plt.xlabel(r'$\Phi=gh$, m$^2$/s$^2$')
    plt.ticklabel_format(axis='both', style='sci')
    plt.ylabel('latitude')
    plt.ticklabel_format(axis='both', style='sci',scilimits=(0,1))
    
    plt.xlim(-10*2, 2*10**5)
    plt.title('t='+str(t*10)+' hours')
    plt.show()
    
    # IMPORTANT ANIMATION CODE HERE
    # Used to keep the limits constant

    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image

Phidata=np.zeros((num_snapshot,J,I))
frms=5

for m in range(num_snapshot): 
    #Phidata[m,:,:]=read_pickle('data/','Phi-'+str(tspaced[m]))
    Phidata[m,:,:]=read_pickle('data/','Phi-'+str(tspaced[m]))
    
write_zonal_gif(mus, Phidata, num_snapshot, frms, 'low-Phibar4-P1-tau10-zonal-Phi.gif')    