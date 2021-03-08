#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:01:55 2021

@author: Home
"""
import numpy as np

def save_output(output,name):
    #store output in correcto format
    data = np.asarray(np.real(output))
    # save to csv file
    np.savetxt(name+'.csv', data, delimiter=',')
    
def load_input(name):
    # import csv file
    data = np.loadtxt(name+'.csv', delimiter=',')
    
    return data
    