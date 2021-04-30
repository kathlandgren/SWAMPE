#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:01:55 2021

@author: Home
"""
import numpy as np

def save_output(output,name):
    #store output in correct format
    data = np.asarray(np.real(output))
    # save to csv file
    np.savetxt(name+'.csv', data, delimiter=',')
    
def load_input(name):
    # import csv file
    data = np.loadtxt(name+'.csv', delimiter=',')
    
    return data
    

def flatten_and_save(output,name):

      
    # reshaping the array from 3D
    # matrice to 2D matrice.
    arr_reshaped = output.reshape(output.shape[0], -1)
      
    # saving reshaped array to file.
    np.savetxt(name+".txt", arr_reshaped)
      
def load_and_restore(name,I):
    
    # retrieving data from file.
    loaded_arr = np.loadtxt(name+'.txt')
      
    # This loadedArr is a 2D array, therefore
    # we need to convert it to the original
    # array shape.reshaping to get original
    # matrice with original shape.
    data = loaded_arr.reshape(
        loaded_arr.shape[0], loaded_arr.shape[1] // I, I)
    
    return data
    