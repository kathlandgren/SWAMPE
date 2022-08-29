# -*- coding: utf-8 -*-
"""
This module contains functions related to generating figures and gifs with SWAMPE.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import os
import imageio


def mean_zonal_wind_plot(plotdata,mus,timestamp,units='hours',customtitle=None,customxlabel=None,savemyfig=False,filename=None,custompath=None,color=None):

    #average data
    zonal_mean=np.mean(plotdata,axis=1)
    
    #convert latitudes to degrees
    Y = np.arcsin(mus)*180/np.pi #generate latitudes in degrees
    
    fig = plt.figure()
    
    if color!=None:
        fig=plt.plot(zonal_mean, Y,color=color)
    else:
        fig=plt.plot(zonal_mean, Y)
    
    if customxlabel==None:    
        plt.xlabel('mean U, m/s')
    else:
        plt.xlabel(customxlabel)
        
    plt.ticklabel_format(axis='both', style='sci')
    plt.ylabel('latitude')
    plt.ticklabel_format(axis='both', style='sci')
    
    if customtitle==None:
        plt.title('Mean zonal winds at '+ str(timestamp) + ' '+units)
    else:
        plt.title(customtitle)
        
    if savemyfig==True:
        if custompath==None:
            path = 'plots/'
            isExist = os.path.exists(path)
            if isExist==False:
                os.mkdir('plots/')
            plt.savefig(path+filename, bbox_inches='tight', dpi=800)
        else:
            plt.savefig(custompath+filename, bbox_inches='tight', dpi=800)
    
    plt.show()
    return fig

#need for colorbar
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def quiver_geopot_plot(U,V,Phi,lambdas,mus,timestamp,sparseness=4,minlevel=None,maxlevel=None,units='hours',customtitle=None,savemyfig=False,filename=None,custompath=None,axlabels=False,colormap=None):
    
    # set up lat-lon grid
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    fig, ax=plt.subplots()
    # Plot the geopotential surface.
    ax.contourf(X, Y, (Phi))

    if minlevel==None:
        minlevel=np.min(Phi)
    if maxlevel==None:
        maxlevel=np.max(Phi)
        
    if colormap==None:
        colormap=cm.nipy_spectral
    #set the colorbar limits
    levels =np.linspace(minlevel, maxlevel)
    CS = plt.contourf(X, Y, (Phi), levels=levels, cmap=colormap, extend='both') 
    colorbar = plt.colorbar(CS,format=ticker.FuncFormatter(fmt))

    # set sparce grid to display wind vector field    
    Xsparse=X[0::sparseness]
    Ysparse=Y[0::sparseness]

    Usparse=U[0::sparseness,0::sparseness]
    Vsparse=V[0::sparseness,0::sparseness]
    
    #to prevent aliasing if saving to pdf
    for c in ax.collections:
        c.set_edgecolor("face")   

    plt.quiver(Xsparse,Ysparse,Usparse,Vsparse)
    
    if axlabels==True:
        ax.set_ylabel('latitude')
        ax.set_xlabel('longitude')
    
    if customtitle==None:
        ax.set_title('Geopotential at '+ str(timestamp) + ' '+units)
    else:
        ax.set_title(customtitle)
    
    if savemyfig==True:
        if custompath==None:
            path = 'plots/'
            isExist = os.path.exists(path)
            if isExist==False:
                os.mkdir('plots/')
            plt.savefig(path+filename, bbox_inches='tight', dpi=800)
        else:
            plt.savefig(custompath+filename, bbox_inches='tight', dpi=800)
        
    return fig
    

def spinup_plot(plotdata,dt,units='hours',customtitle=None,customxlabel=None,customylabel=None,savemyfig=False,filename=None,custompath=None,color=None,legendflag=True,customlegend=None):
    """

    
    """
    tmax=np.shape(plotdata)[0]
    
    if units=='hours':
        tlim=dt*tmax/3600
    elif units=='minutes':
        tlim=dt*tmax/60
    elif units=='seconds':
        tlim=dt*tmax
    else:
        print('Cannot parse units. Acceptable units are: hours, minutes, seconds.')
    
    t = np.linspace(0, tlim, num=tmax, endpoint=True)
    fig_spinup = plt.figure()
    if color!=None:
        fig_spinup=plt.plot(t, plotdata[:,0],color=color[0])
        plt.plot(t, plotdata[:,1],color=color[1])
    else:
        fig_spinup=plt.plot(t, plotdata[:,0])
        plt.plot(t, plotdata[:,1])
    
    
    if legendflag==True:
        if customlegend==None:
            plt.legend(['min winds','RMS winds'],loc="center right")
        else:
            plt.legend(customlegend,loc="center right")
    
    if customxlabel==None:
        plt.xlabel('time, '+units)
    else:
        plt.xlabel(customxlabel)
    
    if customylabel==None:
        plt.ylabel('Winds, m/s')
    else:
        plt.ylabel(customylabel)
        
    plt.ticklabel_format(axis='both', style='sci')
    
    if customtitle==None:
        plt.title('RMS winds during spin-up')
    else:
        plt.title(customtitle)
    
    if savemyfig==True:
        if custompath==None:
            path = 'plots/'
            isExist = os.path.exists(path)
            if isExist==False:
                os.mkdir('plots/')
            plt.savefig(path+filename, bbox_inches='tight', dpi=800)
        else:
            plt.savefig(custompath+filename, bbox_inches='tight', dpi=800)
    
    plt.show()
    return fig_spinup



def gif_helper(fig,dpi=200):
    fig.set_dpi(dpi)
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image

def write_quiver_gif(lambdas,mus,Phidata,Udata,Vdata,timestamps,filename,frms=5,sparseness=4,dpi=200,minlevel=None,maxlevel=None,units='hours',customtitle=None,custompath=None,axlabels=False,colormap=None):
    if minlevel==None:
        minlevel=np.min(Phidata)
    if maxlevel==None:
        maxlevel=np.max(Phidata)
    imageio.mimsave('./'+filename, [gif_helper(quiver_geopot_plot(Udata[i,:,:],Vdata[i,:,:],Phidata[i,:,:],lambdas,mus,timestamps[i],sparseness=sparseness,minlevel=minlevel,maxlevel=maxlevel,units=units,customtitle=customtitle,custompath=None,axlabels=axlabels,colormap=colormap),dpi=200) for i in range(len(timestamps))], fps=frms)
