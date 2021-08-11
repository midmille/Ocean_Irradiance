# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 22:36:14 2020

@author: Miles Miller
"""

"""
This file contains a simple function that solves for the r,w grid of a ROMS output 

"""

from netCDF4 import Dataset 
import os
os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share" 
import seapy 


def make_roms_grid_from_netcdf(file):
    """
    

    Parameters
    ----------
    file : string 
        The netcdf file output from roms. 

    Returns
    -------
    z_r : Array
        The rho grid. The cell centers 
    z_w : Array 
        The w grid. The cell edges. 

    """
        
    roms_data = Dataset(file, 'r')
    
    h = roms_data.variables['h'][:]
    hc = roms_data.variables['hc'][:]
    s_rho = roms_data.variables['s_rho'][:]
    s_w = roms_data.variables['s_w'][:]
    Cs_r = roms_data.variables['Cs_r'][:]
    Cs_w = roms_data.variables['Cs_w'][:]
    zeta = roms_data.variables['zeta'][0,:]
    vtransform = roms_data.variables['Vtransform'][:]
    
    z_r = seapy.roms.depth(vtransform, h, hc, s_rho, Cs_r,zeta)[:,:,:] ##the rho grid 
    z_w = seapy.roms.depth(vtransform, h, hc, s_w, Cs_w,zeta)[:,:,:]  ##the w-grid at eddges 
    
    return z_r, z_w






