
"""
This module contains the functions needed to save and read SWAMP-E data, as well as for continuation.
"""
import pickle
import os

def write_pickle(filename,data,custompath=None):
    """
    
    :param filename: name of the pickle file to be saved
    :type filename: string
    :param data: a Python array of data to be saved
    :type data: any Python type
    :param custompath: path to the custom directory, defaults to None. If None, files will be saved in the parent_directory/data/.
    :type custompath: string, optional

    """
    
    if custompath==None:
        outfile = open('data/'+filename,'wb')
    else:
        outfile = open(custompath+filename,'wb')
        
    pickle.dump(data,outfile)
    outfile.close()
    
def read_pickle(filename,custompath=None):
    """
    
    :param filename: name of the pickle file to be read
    :type filename: string
    :param custompath: path to the custom directory, defaults to None
    :type custompath: string, optional
    :return: any Python type from the pickle file
    :rtype: any Python type

    """
    
    if custompath==None:
        infile = open('data/'+filename,'rb')
    else:
        infile = open(custompath+filename,'rb')
    var = pickle.load(infile)
    infile.close()
    return var

def compute_timestamp(units,t,dt):
    """
    Computes timestamp in appropriate units to append to the saved data files.

    :param units: Units of timestamps on the savefile: 'hours','minutes', or 'seconds'
    :type units: str
    :param t: number of current timestep
    :type t: int
    :param dt: timestep length, in seconds
    :type dt: float
    :return: timestamp in desired units
    :rtype: string

    """
    
    
    if units=='hours':
        timestamp=str(int(dt*t/3600))
    elif units=='minutes':
        timestamp=str(int(dt*t/60))
    elif units=='seconds':
        timestamp=str(int(dt*t))
    else:
        print('Cannot parse units. Acceptable units are: hours, minutes, seconds.')
    
    return timestamp

def compute_t_from_timestamp(units,timestamp,dt):
    """_summary_

    :param units: Units of timestamps on the savefile: 'hours','minutes', or 'seconds'
    :type units: str
    :param timestamp: Timestamp in specified units
    :type timestamp: int
    :param dt: timestep length, in second s
    :type dt: float

    return: number of timestep to continue the simulation
    rtype: int
    """

    if units=='hours':
        t=int(timestamp*3600/dt) 
    elif units=='minutes':
        t=int(timestamp*60/dt)
    elif units=='seconds':
        t=int(timestamp/dt)
    else:
        print('Cannot parse units. Acceptable units are: hours, minutes, seconds.')

    return t 

def save_data(timestamp,etadata, deltadata, Phidata, U,V, spinupdata,geopotdata,  custompath=None):
    """
    Saves the data for plotting and continuation purposes.
    
    :param timestamp: timestamp to be used for naming saved files
    :type timestamp: string
    :param etadata: python array of the data for absolute vorticity eta
    :type etadata: array of float JxI
    :param deltadata: python array of the data for divergence delta
    :type deltadata: array of float JxI
    :param Phidata: python array of the data for geopotential Phi
    :type Phidata: array of float JxI
    :param U: python array of the data for zonal winds U
    :type U: array of float JxI
    :param V: python array of the data for meridional winds V
    :type V: array of float JxI
    :param spinupdata: time series array of minimum length of wind vector and RMS winds
    :type spinupdata: array of float
    :param geopotdata: time series array of minimum and maximum of the geopotential Phi
    :type geopotdata: array of float
    :param custompath: path to the custom directory, defaults to None. If None, files will be saved in the parent_directory/data/
    :type custompath: string, optional


    """
    
    if custompath==None:
        path = 'data/'
        isExist = os.path.exists(path)
        if isExist==False:
            os.mkdir('data/')

        write_pickle('eta-'+timestamp,  etadata) 
        write_pickle('delta-'+timestamp,  deltadata) 
        write_pickle('Phi-'+timestamp,  Phidata)   
        write_pickle('U-'+timestamp, U)    
        write_pickle('V-'+timestamp, V) 
        
        write_pickle('spinup-winds', spinupdata) 
        write_pickle('spinup-geopot', geopotdata) 
    
    else:
        
        write_pickle('eta-'+timestamp,  etadata,custompath=custompath) 
        write_pickle('delta-'+timestamp,  deltadata,custompath=custompath) 
        write_pickle('Phi-'+timestamp,  Phidata, custompath=custompath)   
        write_pickle('U-'+timestamp, U, custompath=custompath)    
        write_pickle('V-'+timestamp, V, custompath=custompath) 
        
        write_pickle('spinup-winds', spinupdata, custompath=custompath) 
        write_pickle('spinup-geopot', geopotdata, custompath=custompath) 
        
def load_data(timestamp,custompath=None):
    """
    
    :param timestamp: timestamp used for naming saved files
    :type timestamp: string
    :param custompath: path to the custom directory, defaults to None. If None, files will be saved in the parent_directory/data/
    :type custompath: string
    
    :return: arrays of eta, delta, Phi, U, V
    :rtype: arrays of float JxI

    """
    
    if custompath==None:
        eta=read_pickle('eta-'+str(timestamp))
        delta=read_pickle('delta-'+str(timestamp))
        Phi=read_pickle('Phi-'+str(timestamp))
        U=read_pickle('U-'+str(timestamp))
        V=read_pickle('V-'+str(timestamp))
    else:
        eta=read_pickle('eta-'+str(timestamp),custompath=custompath)
        delta=read_pickle('delta-'+str(timestamp),custompath=custompath)
        Phi=read_pickle('Phi-'+str(timestamp),custompath=custompath)
        U=read_pickle('U-'+str(timestamp),custompath=custompath)
        V=read_pickle('V-'+str(timestamp),custompath=custompath)
        
    return eta, delta, Phi, U, V 
    