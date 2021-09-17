# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 11:49:27 2021

@author: miles

Read and Play with Ocean Color data. 
"""

import numpy as np


def Make_Data_Dict(file, field_names, skip_header):
    """
    This reads in a CSV data file from the ocean color website and outputs a
    dictionary with keys ccorrespodinbg to the name of each field and values corresponding 
    to the data array of each field. 

    Parameters
    ----------
    file : String
        Path to CSV data file.
    field_names : String
        The variable names as a comma seperated string. This should be copied from the 
        csv file.
    skip_header : Integer
        The numnber of header lines to skip in the file.

    Returns
    -------
    data_dict : Dictionary
        The data dictionary output. Keys are the field names. Values are the data 
        arrays. The data type of the value arrays are of np.dtype = 'object'.

    """
    

    ## This makes a size 2278 numpy array. 
    ## Each index has a length 22 tuple  which contains each field.  
    data = np.genfromtxt(file, dtype=None, skip_header = skip_header, delimiter=',')
    
    fields = field_names.split(',')
    
    data_dict = {}
    
    ## Creating the empty arrays for each field in dictionary.
    for j,field in enumerate(fields):
        data_dict[field] = np.empty(len(data), dtype='object')
        # if type(data[0][j] == float) or type(data[0][j] == int): 
        #     data_dict[field] = np.empty(len(data), dtype='s10')
        # else:
        #     data_dict[field] = np.empty(len(data), dtype='s10')
    
    ## Loop over each tuple/row in data. 
    for k,row in enumerate(data): 
        ## Checking for compatible sizes
        assert len(row) == len(fields)
        ## Loop over each field and store corresponidng data. 
        for j,field in enumerate(fields):
            ## for each dictionary value, the correct index of data field array 
            ## is stored with corresponding data value from current tuple/row of data set.
            data_dict[field][k] = row[j]
            
    return data_dict


    
    

    

    
    
    
    
    
    
    
    
