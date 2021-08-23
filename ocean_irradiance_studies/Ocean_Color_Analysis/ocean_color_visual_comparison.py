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
import numpy as np



def Plot_Comparison(chl_mat, chl_irr, chl_obs, obs_Ygrid, obs_Xgrid):
    
    fig, axes = plt.subplots(1,3)
    (ax1,ax2,ax3) = axes 
    
    vmin,vmax = (0,1)
    cmap='nipy_spectral'
    
    im1 = ax1.pcolormesh(chl_mat, vmin=vmin, vmax=vmax, cmap=cmap)
    im2 = ax2.pcolormesh(chl_irr, vmin=vmin, vmax=vmax, cmap=cmap)
    
    ## plotting the observed as a scatter plot. 
    s=2
    ax3.scatter(obs_Xgrid, obs_Ygrid, c = chl_obs, s=s, cmap = cmap, vmin=vmin, vmax=vmax)
    ax3.set_ylim([0,181])
    ax3.set_xlim([0,186])
    #ax2.scatter(obs_Xgrid, obs_Ygrid, c = chl_obs, s=s, cmap = cmap, vmin=vmin, vmax=vmax)

    fig.colorbar(im1, ax=axes.ravel().tolist(), label=r'Chlorophyll [$\frac{\mathrm{mg chl_a}}{m^3}}$]' ) 

    ## Labels 
    ax1.set_title('P. Mattern Chlorophyll Model')
    ax2.set_title('Irradiance Chlorophyll Model')
    ax3.set_title('Observed Data')

    for ax in axes: 
        ax.set_xticks([])
        ax.set_yticks([])


    fig.show()
    return 


def Chlorophyll_Concentration_Comparison(irradiance_nc_file, mattern_nc_file, obs_nc_file):
    
    ## Current time step.
    stp = 2
    
    ## Paul's chlorophyll.
    chl_mat = Dataset(mattern_nc_file).variables['chlorophyll'][stp,-1,:,:]
    
    ## Irradiance chlorophyll
    R_nc = RRO.ROMS_netcdf(irradiance_nc_file, Init_ROMS_Irr_Params=(True))
    chl_irr = R_nc.OCx[stp,:,:]
    ## Using the mask to make the land nan for better plotting
    chl_irr[R_nc.maskr == 0] = np.NaN
    
    ## Observed Chlorophyll
    ## The limit of obeserved values to be plotted in grid cell points I think. 
    z_lim = 41

    obs_nc = Dataset(obs_nc_file)
    obs_Zgrid = obs_nc.variables['obs_Zgrid'][:]

    chl_obs = obs_nc.variables['obs_value'][obs_Zgrid > z_lim]
    obs_Ygrid = obs_nc.variables['obs_Ygrid'][obs_Zgrid > z_lim]
    obs_Xgrid = obs_nc.variables['obs_Xgrid'][obs_Zgrid > z_lim]


    
    ## Plotting these different concentrations
    Plot_Comparison(chl_mat, chl_irr, chl_obs, obs_Ygrid, obs_Xgrid)
    
    
    return 


if __name__ == '__main__': 


    import argparse
    
    parser = argparse.ArgumentParser(description='Plotting Some Chlorophyll Comparisons')
    parser.add_argument('irradiance_nc_file', help = "Complete Path to irradiance nc file" )
    parser.add_argument('mattern_nc_file', help = "Complete Path to P. Mattern's model output" )
    parser.add_argument('obs_nc_file', help = "Complete Path to P. Mattern's observed chlorophyll file" )
    # parser.add_argument('dest_file', help='Path to Destination Directory. Saved as pickle')
    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    args = parser.parse_args()
    
    Chlorophyll_Concentration_Comparison(args.irradiance_nc_file, args.mattern_nc_file, args.obs_nc_file)




































