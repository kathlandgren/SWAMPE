
"""
Created on Mon Jun 22 13:57:02 2020

@author: ek672
This script contains the plotting functions
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import imageio
import matplotlib.ticker as ticker
import matplotlib.colors as colors

def spectral_plot(plotdata):
    """Generates plots for the spectral space with dimensions M, N

    :param plotdata: first variable, MxN 
    :type statevar1: float

    
    """
    ## Plot the approximation 
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # Get the dimensions
    size = plotdata.shape
    M = size[0]
    N = size[1]
    
    # Make data.
    X = np.arange(N)
    Y = np.arange(M)
    X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, (np.real(plotdata)), cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    
    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    # Add axis labels
    ax.set_xlabel('N', fontsize=20)
    ax.set_ylabel('M', fontsize=20)
    
    plt.show()
    
def fourier_plot(plotdata,mus):
    """Generates plots for the fourier space with dimensions J, M

    :param plotdata: first variable, JxM+1 
    :type statevar1: float

    
    """
    ## Plot the approximation 
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # Get the dimensions
    size = plotdata.shape
    M = size[1]
    
    # Make data.
    X = np.arange(M)
    Y = mus
    X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, (np.real(plotdata)), cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    
    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    # Add axis labels
    ax.set_xlabel('M', fontsize=20)
    ax.set_ylabel('$\mu$', fontsize=20, rotation = 0)
    
    plt.show()
    
def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)    
    
def physical_plot(plotdata,mus,lambdas):
    """Generates plots for the fourier space with dimensions J, M

    :param plotdata: first variable, JxM+1 
    :type statevar1: float

    
    """
    ## Plot the approximation 
    fig= plt.figure()
    
    # Make data.
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    cp = plt.contourf(X, Y, plotdata,30)
    norm = mpl.colors.Normalize(vmin=0, vmax=5*10**6)
    cb = plt.colorbar(cp,format=ticker.FuncFormatter(fmt),norm=norm)
    #cb.set_clim(0, 5*10**(6))
    
    plt.show()
    
    


def plot_for_offset(t, z_max,lambdas,mus,Xidata):
    
    """
    Generates the plot to generate the gif later
    
    :param t: time, t<tmax
    :type t: int
    :param z_max: max height of the plot
    :type z_max: float
    :param lambdas: longitudes
    :type lambdas: list
    :param mus: Gaussian latitudes
    :type mus: list
    :param Xidata: state variable data for all time
    :type Xidata: list
    """
    # Make data.
    X = lambdas
    Y = mus
    X, Y = np.meshgrid(X, Y)

    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Plot the surface.
    surf = ax.plot_surface(X, Y, np.real(Xidata[t,:,:]), cmap=cm.coolwarm, linewidth=0, antialiased=False)

    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    #fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
    # IMPORTANT ANIMATION CODE HERE
    # Used to keep the limits constant
    ax.set_zlim(-z_max, z_max)

    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image



def write_gif( z_max,lambdas,mus,Xidata,tmax,frms,string):
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
    imageio.mimsave('./'+string, [plot_for_offset(i, z_max,lambdas,mus,Xidata) for i in range(tmax)], fps=frms)
    

def plot_for_offset_2D(t,lambdas,mus,plotdata):
    
    """
    Generates the plot to generate the gif later
    
    :param t: time, t<tmax
    :type t: int
    :param z_max: max height of the plot
    :type z_max: float
    :param lambdas: longitudes
    :type lambdas: list
    :param mus: Gaussian latitudes
    :type mus: list
    :param Xidata: state variable data for all time
    :type Xidata: list
    """
    # Make data.
       ## Plot the approximation 
    fig,ax= plt.subplots(figsize=(10,5))
    
    # Make data.
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    cp = plt.contourf(X, Y, plotdata[t,:,:],30)

    cb = plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
    #cb.set_clim(0, 5*10**(6))
    
    plt.title('t='+str(t))
    plt.show()
    # IMPORTANT ANIMATION CODE HERE
    # Used to keep the limits constant

    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image

def write_2D_gif(lambdas,mus,Xidata,tmax,frms,string):
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
    imageio.mimsave('./'+string, [plot_for_offset_2D(i,lambdas,mus,Xidata) for i in range(tmax)], fps=frms)
    
def quiver_for_offset_2D(U,V,Phi,lambdas,mus,t,sparseness,test,a1):
    U=U[t,:,:]
    V=V[t,:,:]
    Phi=Phi[t,:,:]
    
    
    fig,ax= plt.subplots(figsize=(10,5))
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    #X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    plt.contourf(X, Y, Phi,30)
    cb = plt.colorbar(format=ticker.FuncFormatter(fmt))
    
    Xsparse=X[0::sparseness]
    Ysparse=Y[0::sparseness]
    #norm=np.sqrt(np.multiply(U,U)+np.multiply(V,V))
    Usparse=U[0::sparseness,0::sparseness]
    Vsparse=V[0::sparseness,0::sparseness]
    # plt.quiver(x, y[skip], u[skip], v[skip], color='black', headwidth=1, scale = 10, headlength=4)
    
    # ax.quiver(X,Y,U,V)
    # Xsparse, Ysparse = np.meshgrid(Xsparse, Ysparse)
    plt.quiver(Xsparse,Ysparse,Usparse,Vsparse)
    plt.title('t='+str(t)+', test='+str(test)+', alpha='+str(a1))
    plt.show()
    # IMPORTANT ANIMATION CODE HERE
    # Used to keep the limits constant

    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image


def write_quiver_gif(lambdas,mus,Phidata,Udata,Vdata,tmax,frms,string,sparseness,test,a1,minlevel,maxlevel):
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
    imageio.mimsave('./'+string, [quiver_for_offset_2D(Udata,Vdata,Phidata,lambdas,mus,i,sparseness,test,a1,minlevel,maxlevel) for i in range(tmax)], fps=frms)
# # def gif_contour():
# #     # Generate grid for plotting
# #     X = lambdas
# #     Y = mus
# #     X, Y = np.meshgrid(X, Y)
    
# #     fig = plt.figure()
# #     plt.xlabel(r'longitude')
# #     plt.ylabel(r'latitude')

# #     # animation function
# #     z = var[i,:,0,:].T
# #     cont = plt.contourf(x, y, z, 25)
# #         if (tslice == 0):
# #             plt.title(r't = %1.2e' % t[i] )
# #         else:
# #             plt.title(r't = %i' % i)
    
# #         return cont  
    
# #     anim = animation.FuncAnimation(fig, animate, frames=Nt)
    
# #     anim.save('animation.mp4')
    
def quiver_plot(U,V,lambdas,mus,sparseness):
    
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    Xsparse=X[0::sparseness]
    Ysparse=Y[0::sparseness]
    #norm=np.sqrt(np.multiply(U,U)+np.multiply(V,V))
    Usparse=U[0::sparseness,0::sparseness]
    Vsparse=V[0::sparseness,0::sparseness]
    # plt.quiver(x, y[skip], u[skip], v[skip], color='black', headwidth=1, scale = 10, headlength=4)
    
    fig, ax = plt.subplots()
    # ax.quiver(X,Y,U,V)
    ax.quiver(Xsparse,Ysparse,Usparse,Vsparse)
    
    plt.show()
    

    
def quiver_geopot_plot(U,V,Phi,lambdas,mus,t,sparseness,test,a1,minlevel,maxlevel):
    
    X = lambdas*180/np.pi
    Y = np.arcsin(mus)*180/np.pi
    #X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.

    plt.contourf(X, Y, (Phi))
    #plt.colorbar(extend='both')

    levels =np.linspace(minlevel, maxlevel) #set the colorbar limits
    CS = plt.contourf(X, Y, np.log10(Phi), levels=levels, cmap=cm.jet, extend='both')
    
    colorbar = plt.colorbar(CS)

    #cb = plt.colorbar(format=ticker.FuncFormatter(fmt),extend='both')
 

    #plt.colorbar(extend='both')
    #plt.clim(0, 10**4)
    
    Xsparse=X[0::sparseness]
    Ysparse=Y[0::sparseness]
    #norm=np.sqrt(np.multiply(U,U)+np.multiply(V,V))
    Usparse=U[0::sparseness,0::sparseness]
    Vsparse=V[0::sparseness,0::sparseness]
    # plt.quiver(x, y[skip], u[skip], v[skip], color='black', headwidth=1, scale = 10, headlength=4)
    
    # ax.quiver(X,Y,U,V)
    # Xsparse, Ysparse = np.meshgrid(Xsparse, Ysparse)
    plt.quiver(Xsparse,Ysparse,Usparse,Vsparse)
    plt.title('t='+str(t)+', test='+str(test)+', alpha='+str(a1))
    plt.show()
        
#     fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.plot(series1_x, series1_y)

# ax2 = fig.add_subplot(111)
# ax2.plot(series2_x, series2_y)

# ax3 = fig.add_subplot(111)
# ax3.scatter(series2_x, series2_y)
    