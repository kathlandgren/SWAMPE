
"""
This module contains the files needed to 
"""
import pickle


def write_pickle(filename,data):
    outfile = open('data/'+filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()
    
def read_pickle(filename):
    infile = open('data/'+filename,'rb')
    var = pickle.load(infile)
    infile.close()
    return var