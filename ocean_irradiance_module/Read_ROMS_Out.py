# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:27:43 2021

@author: miles
"""

## User created modules
# from Ocean_Irradiance_ROMS import make_roms_grid_from_netcdf 
## Python modules
from netCDF4 import Dataset
from ocean_irradiance_module.PARAMS import Param_Init 




class ROMS_netcdf:
    """
    
    """
    
    def __init__(self, file, Init_ROMS_Irr_Params=False):
        
        self.file = file
        PI = Param_Init(self.file)
        ##NOTE: this takes a little while to compute, save object as pickle for recurring problem. 
        ## The ROMS data as a data set object
        roms_nc = Dataset(self.file, 'r')
        
        ## The ROMS grid as defined by coordinates at rho or velcoties(w).
        # self.z_r,self.z_w = make_roms_grid_from_netcdf(self.file) 
        ## Time in seconds since 1999, 01/01, 00:00:003.
        self.ocean_time = roms_nc.variables['ocean_time'][:]
        ## The shape of the phytoplankton concentration arrays.
        self.nti,self.nzi,self.nyi,self.nxi = roms_nc.variables["diatom"].shape 
        ## Getting diatom array.
        self.diatom =  roms_nc.variables["diatom"][:]
        ## Getting nanophytoplankton array.
        self.nanophyt =  roms_nc.variables["nanophytoplankton"][:]
        ## Changing to chl
        self.chl_diatom = PI.Chl2NL * self.diatom
        ## Changing to chl
        self.chl_nanophyt = PI.Chl2NS * self.nanophyt
        ## The mask that defines land or not. Land is taken as 0, while water = 1. 
        self.maskr = roms_nc.variables['mask_rho'][:] 
        
        ## ROMS grids 
        self.z_r = roms_nc.variables['z_rho'][:]
        self.z_w = roms_nc.variables['z_w'][:]      
        ## Irradiance Arrays
        self.Ed1 = roms_nc.variables['Ed_1'][:]
        self.Ed2 = roms_nc.variables['Ed_2'][:]
        self.Es1 = roms_nc.variables['Es_1'][:]
        self.Es2 = roms_nc.variables['Es_2'][:]
        self.Eu1 = roms_nc.variables['Eu_1'][:]
        self.Eu2 = roms_nc.variables['Eu_2'][:]
        ## Irradiance Grid
        self.z_irr1 = roms_nc.variables['z_irr_1'][:]
        self.z_irr2 = roms_nc.variables['z_irr_2'][:]
        ## Ocean Color output
        self.OCx = roms_nc.variables['OCx'][:]
        
        ## The date of the given ROMS output. 
        # self.roms_date = seapy.roms.num2date(nc = self.roms_data)
        
        if Init_ROMS_Irr_Params: 
    
            ## This reads in irradiance params from ROMS nc out only if there are any...
            ## otherwise its best to use PARAMS class. 
            self.Ed0 = roms_nc.variables['Ed0'][:]
            self.Es0 = roms_nc.variables['Es0'][:]
            self.Euh = roms_nc.variables['Euh'][:]
            self.wavelengths = roms_nc.variables['wavelengths'][:]
            self.N_irr = roms_nc.variables['N_irr'][:]
            
            ## Getting coefficients and putting them into a form more readble by python verion of irradiance
            a_wat_lam = roms_nc.variables['a_wat_lam']
            b_wat_lam = roms_nc.variables['b_wat_lam']
            a_diatom_lam = roms_nc.variables['a_diatom_lam']
            b_diatom_lam = roms_nc.variables['b_diatom_lam']
            a_nanophyt_lam = roms_nc.variables['a_nanophyt_lam']
            b_nanophyt_lam = roms_nc.variables['b_nanophyt_lam']
            
            self.ab_wat = {}
            self.ab_diat = {}
            self.ab_syn = {} 
            
            for k, lam in enumerate(self.wavelengths):
                self.ab_wat[lam] = (a_wat_lam[k], b_wat_lam[k])
                self.ab_diat[lam] = (a_diatom_lam[k], b_diatom_lam[k])
                self.ab_syn[lam] = (a_nanophyt_lam[k], b_nanophyt_lam[k])
                










##OLD STUFF 
##--------------
# import os
## Weird dictionary edit necessary because my basemap module lacks the following key:value pair. 
# os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share" 

## Having trouble importing seapy at the moment. 
# import seapy 


# def make_roms_grid_from_netcdf(file):
#     """
#     Retreives the ROMS grid from file. 
    

#     Parameters
#     ----------
#     file : string 
#         The netcdf file output from roms. 

#     Returns
#     -------
#     z_r : Array
#         The rho grid. The cell centers 
#     z_w : Array 
#         The w grid. The cell edges. 

#     """
    
#     roms_data = Dataset(file, 'r')
    
#     h = roms_data.variables['h'][:]
#     hc = roms_data.variables['hc'][:]
#     s_rho = roms_data.variables['s_rho'][:]
#     s_w = roms_data.variables['s_w'][:]
#     Cs_r = roms_data.variables['Cs_r'][:]
#     Cs_w = roms_data.variables['Cs_w'][:]
#     zeta = roms_data.variables['zeta'][0,:]
#     vtransform = roms_data.variables['Vtransform'][:]
    
#     z_r = seapy.roms.depth(vtransform, h, hc, s_rho, Cs_r,zeta)[:,:,:] ##the rho grid 
#     z_w = seapy.roms.depth(vtransform, h, hc, s_w, Cs_w,zeta)[:,:,:]  ##the w-grid at eddges 

#     return z_r, z_w

        
   
