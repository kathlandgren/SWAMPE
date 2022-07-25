
"""
This module contains the files needed to 
"""
import pickle
import os

def write_pickle(filename,data,custompath=None):
    
    if custompath==None:
        outfile = open('data/'+filename,'wb')
    else:
        outfile = open(custompath+filename,'wb')
        
    pickle.dump(data,outfile)
    outfile.close()
    
def read_pickle(filename,custompath=None):
    
    if custompath==None:
        infile = open('data/'+filename,'rb')
    else:
        infile = open(custompath+filename,'rb')
    var = pickle.load(infile)
    infile.close()
    return var

def save_data(timestamp, Phidata, U,V, spinupdata,geopotdata,  custompath=None):
    
    print(custompath)
    
    if custompath==None:
        path = 'data/'
        isExist = os.path.exists(path)
        if isExist==False:
            os.mkdir('data/')

        write_pickle('Phi-'+timestamp,  Phidata)   
        write_pickle('U-'+timestamp, U)    
        write_pickle('V-'+timestamp, V) 
        
        write_pickle('spinup-winds', spinupdata) 
        write_pickle('spinup-geopot', geopotdata) 
    
    else:
        
        write_pickle('Phi-'+timestamp,  Phidata, custompath=custompath)   
        write_pickle('U-'+timestamp, U, custompath=custompath)    
        write_pickle('V-'+timestamp, V, custompath=custompath) 
        
        write_pickle('spinup-winds', spinupdata, custompath=custompath) 
        write_pickle('spinup-geopot', geopotdata, custompath=custompath) 