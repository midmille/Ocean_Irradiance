# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 17:10:18 2021

@author: Miles Miller

This file defines all variables and main constants. 
"""

## User created modules
import ROMS_grid

## Python modules
from netCDF4 import Dataset 
import os
## Weird dictionary edit necessary because my basemap module lacks the following key:value pair. 
os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share" 
import seapy 

class Param_Init:
    """
    This class defines the default constants and ROMS parameters.
    """
    def __init__(self, file):
        #####################################INITIAL PARAMS###########################
        ##############################################################################
        # file = '/Users/Miles Miller/da_fwd_002.nc'
        ## The main ROMS output netCDF file is taken as an initial input of class.
        self.file = file
        ## The default boundary conditions of problem following John Wilkin.
        self.E_d_0 = .7
        self.E_s_0 = 1 - self.E_d_0
        self.E_u_h = 0 
        
        ## Wavelengths taken from John Wilkin email regarding satelite wavelengths.
        self.wavelengths = [410,443,486,551,638,671] 
        ## The ROMS data as a data set object
        self.roms_data = Dataset(self.file, 'r')
        ## The Chlorophyll to Nitrogen unit change. 
        self.Chl2NL = 1.59
        self.Chl2NS = .7950
        ## The time step of the problem
        self.time_step_index = 0 
        
        ##NOTE: this takes a little while to compute, save object as pickle for recurring problem. 
        
        ## The ROMS grid as defined by coordinates at rho or velcoties(w).
        self.z_r,self.z_w = ROMS_grid.make_roms_grid_from_netcdf(self.file) 
        ## Time in seconds since 1999, 01/01, 00:00:003.
        self.ocean_time = self.roms_data.variables['ocean_time'][:]
        ## The shape of the phytoplankton concentration arrays.
        self.nti,self.nzi,self.nyi,self.nxi = self.roms_data.variables["diatom"].shape 
        ## Getting diatom array.
        self.diatom =  self.roms_data.variables["diatom"][self.time_step_index,:,:,:]
        ## Getting nanophytoplankton array.
        self.nanophyt =  self.roms_data.variables["nanophytoplankton"][self.time_step_index,:,:,:]
        ## Changing to chl
        self.chl_diatom = self.Chl2NL * self.diatom
        ## Changing to chl
        self.chl_nanophyt = self.Chl2NS * self.nanophyt
        ## The mask that defines land or not. Land is taken as 0, while water = 1. 
        self.maskr = self.roms_data.variables['mask_rho'][:] 
        ## The date of the given ROMS output. 
        self.roms_date = seapy.roms.num2date(nc = self.roms_data)
        
        
        ##############################################################################
        ##############################################################################
        
    
