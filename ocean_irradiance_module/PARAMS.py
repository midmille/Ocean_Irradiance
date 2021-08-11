# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 17:10:18 2021

@author: Miles Miller

This file defines all variables and main constants. 
"""






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
        ## The Chlorophyll to Nitrogen unit change. 
        self.Chl2NL = 1.59
        self.Chl2NS = .7950

        ##############################################################################
        ##############################################################################
        
    
