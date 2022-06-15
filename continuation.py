#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:01:55 2021

@author: Home
"""
import pickle

# def save_output(output,name):
#     """
#     Saves the output as csv

#     Parameters
#     ----------
#     :param output: 1D or 2D array to save
#     :type output: array of float64
    
#     :param name: file name
#     :type name: string

#     """
#     #store output in correct format
#     data = np.asarray(np.real(output))
#     # save to csv file
#     np.savetxt(name+'.csv', data, delimiter=',')
    
# def load_input(name):
#     """
#     Loads data from a .csv file

#     Parameters
#     ----------
#     :param name: file name 
#     :type name: string

#     Returns
#     -------
#     :return data: data from file
#     :rtype: 1D or 2D array of float64
#     """
#     # import csv file
#     data = np.loadtxt(name+'.csv', delimiter=',')
    
#     return data
    

# def flatten_and_save(output,name):
#     """
#     Reshapes a 3D array into 2D and saves as .txt
    
#     Paramerers
#     ----------
#     :param output: 3D array
#     :type output: array of float64
    
#     :param name: file name
#     :type name: string
    
#     """
      
#     # reshaping the array from 3D
#     # matrice to 2D matrice.
#     arr_reshaped = output.reshape(output.shape[0], -1)
      
#     # saving reshaped array to file.
#     np.savetxt(name+".txt", arr_reshaped)
      
# def load_and_restore(name,I):
#     """
#     Imports 2D array and reshapes it back into 3D.
    
#     Paramerers
#     ----------

#     :param name: file name
#     :type name: string
    
#     :param I: number of longitudes
#     :type I: int
    
#     Returns
#     -------
    
#     :return data: 3D data array
#     :rtype: array of float64
    
#     """
          
    
#     # retrieving data from file.
#     loaded_arr = np.loadtxt(name+'.txt')
      
#     # This loadedArr is a 2D array, therefore
#     # we need to convert it to the original
#     # array shape.reshaping to get original
#     # matrice with original shape.
#     data = loaded_arr.reshape(
#         loaded_arr.shape[0], loaded_arr.shape[1] // I, I)
    
#     return data
    


def write_pickle(filename,data):
    outfile = open('data/'+filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()
    
def read_pickle(filename):
    infile = open('data/'+filename,'rb')
    var = pickle.load(infile)
    infile.close()
    return var