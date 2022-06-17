#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:01:55 2021

@author: Home
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