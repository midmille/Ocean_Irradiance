# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:17:11 2021

@author: miles

This file creates a visual comparison of my ocean color chlorphyll model 
to Paul Mattwerns model and to data. 

"""
from netCDF4 import Dataset
import ocean_irradiance_module.Read_ROMS_Out as RRO 
import matplotlib.pyplot as plt 



def Plot_Comparison(chl_mat, chl_irr, chl_obs, obs_Ygrid, obs_Xgrid):
    
    fig, (ax1,ax2) = plt.subplots(1,2)
    
    vmin,vmax = (0,20)
    cmap='hsv'
    
    im1 = ax1.pcolormesh(chl_mat, vmin=vmin, vmax=vmax, cmap=cmap)
    im2 = ax2.pcolormesh(chl_irr, vmin=vmin, vmax=vmax, cmap=cmap)
    
    ## plotting the observed as a scatter plot. 
    ax1.scatter(obs_Xgrid, obs_Ygrid, c = chl_obs, cmap = cmap)

    return 


def chlorphyll_concentration_comparison(irradiance_nc_file, mattern_nc_file, obs_nc_file):
    
    ## Current time step.
    stp = 2
    
    ## Paul's chlorophyll.
    chl_mat = Dataset(mattern_nc_file).variables['chlorophyll'][stp,-1,:,:]
    
    ## Irradiance chlorophyll
    R_nc = RRO.ROMS_netcdf(irradiance_nc_file, Init_ROMS_Irr_Params=(True))
    chl_irr = R_nc.OCx[stp,:,:]
    
    ## Observed Chlorophyll
    obs_nc = Dataset(obs_nc_file)
    chl_obs = obs_nc.variables['obs_value'][:]
    obs_Ygrid = obs_nc.variables['obs_Ygrid'][:]
    obs_Xgrid = obs_nc.variables['obs_Xgrid'][:]
    
    ## Plotting these different concentrations
    Plot_Comparison()
    
    
    return 
    