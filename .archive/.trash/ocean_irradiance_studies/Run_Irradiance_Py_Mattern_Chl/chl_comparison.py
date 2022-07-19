# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 08:47:14 2021

@author: miles

ABOUT 
-----

"""

## Global mods
import numpy as np
import matplotlib.pyplot as plt 
import os
import pickle
import sys
from netCDF4 import Dataset

## appending ocean _irradiance module to path
ocean_irradiance_module_path = os.path.abspath('../..')
sys.path.append(ocean_irradiance_module_path)

##User mods
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
from ocean_irradiance_module.Ocean_Irradiance_ROMS import Ocean_Irradiance_Field 
from ocean_irradiance_module.Ocean_Irradiance_ROMS import ocean_color_sol



def Irradiance_Field_py(R_nc, M_nc, nstp):
    """
    

    Parameters
    ----------
    R_nc : ROMS_netcdf object
        Holds all the returned arrays from ROMS output and paramters used in the simulation. 
    M_nc : Mattern Netcdf4 Dataset object.
    nstp : Float
        The given time step to calculate for.

    Returns
    -------
    irr_py : Dictionary
        The keys are the wavelengths and the values are the Eu at surface arrays.
        This is calculated by python code.
    irr_ROMS : Dictionary
        The keys are the wavelengths and the values are the Eu at surface arrays. 
        This is calculated from within the ROMS Fortran code.
    """
    
    ## The name of the file that the python Eu dict will be saved to as pickle.
    save_file = f'irr_dict_nstp_{nstp}.p'
    save_dir = 'Irr_Field_Out'
    save_path = f'{save_dir}/{save_file}'
    ## Python calculated Eu at surface 
    ##---------------------------------
    ## Checking if save file exists
    ## if it doesn't exist then redo calculation
    if os.path.exists(save_path) == False:
        ## User input required so not to overwrite or redo unnecessarily... 
        y_or_n = input('File does not exist, continue with calculation? [y/n] ')
        if y_or_n == 'n': 
            sys.exit('Stopping...')
        elif y_or_n == 'y':
            print('ok, statrting calculations...')
            
        mask = np.ones((R_nc.nyi, R_nc.nxi))
        irr_field_py = {}
        for lam in R_nc.wavelengths:
            print('Current Wavelength:', lam)
            irr_field_py[lam] = Ocean_Irradiance_Field(mask, 
                                              R_nc.ab_wat[lam], 
                                              R_nc.ab_diat[lam], 
                                              R_nc.ab_syn[lam], 
                                              M_nc['chlL'][nstp,:,:,:], 
                                              M_nc['chlS'][nstp,:,:,:], 
                                              R_nc.z_r[nstp,:,:,:], 
                                              R_nc.Ed0, 
                                              R_nc.Es0, 
                                              R_nc.Euh,
                                              N= R_nc.N_irr)
        
        pickle.dump(irr_field_py, open(save_path, "wb"))
        print('Python calculation complete and saved')
        
    ## if the save file does exist then just load it. gity 
    elif os.path.exists(save_path) == True:
        print(f'Irradiance save file exists! Loading python calculated irradiance field from file "{save_file}"...')
        irr_field_py = pickle.load(open(save_path,'rb'))
        print('Yay, file loaded :)')
        

    return irr_field_py



def Plot_Comparison(chl_mat, chl_irr, chl_py, chl_obs, obs_Ygrid, obs_Xgrid):
    
    fig, axes = plt.subplots(2,2)
    ax1 = axes[0,0]
    ax2 = axes[0,1]
    ax3 = axes[1,0]
    ax4 = axes[1,1] 

    vmin,vmax = (0,1)
    cmap='nipy_spectral'
    
    im1 = ax1.pcolormesh(chl_mat, vmin=vmin, vmax=vmax, cmap=cmap)
    im2 = ax2.pcolormesh(chl_irr, vmin=vmin, vmax=vmax, cmap=cmap)
    im3 = ax3.pcolormesh(chl_py, vmin=vmin, vmax=vmax, cmap=cmap)
    
    
    ## plotting the observed as a scatter plot. 
    s=2
    ax4.scatter(obs_Xgrid, obs_Ygrid, c = chl_obs, s=s, cmap = cmap, vmin=vmin, vmax=vmax)
    ax4.set_ylim([0,181])
    ax4.set_xlim([0,186])
    #ax2.scatter(obs_Xgrid, obs_Ygrid, c = chl_obs, s=s, cmap = cmap, vmin=vmin, vmax=vmax)

    fig.colorbar(im1, ax=axes.ravel().tolist(), label=r'Chlorophyll [$\frac{\mathrm{mg chl_a}}{m^3}}$]' ) 

    ## Labels 
    ax1.set_title('P. Mattern Chlorophyll Model')
    ax2.set_title('Irradiance Chlorophyll Model')
    ax3.set_title('Python Irradiance with Mattern Chl')
    ax4.set_title('Observed Data')

    for ax in axes.ravel().tolist(): 
        ax.set_xticks([])
        ax.set_yticks([])


    fig.show()
    return 


def Chlorophyll_Concentration_Comparison(stp, R_nc, M_nc, obs_nc_file, chl_py):
    
    ## Paul's chlorophyll.
    chl_mat = M_nc.variables['chlorophyll'][stp,-1,:,:]
    
    ## Irradiance chlorophyll
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
    Plot_Comparison(chl_mat, chl_irr, chl_py, chl_obs, obs_Ygrid, obs_Xgrid)
    
    
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

    ## The ROMS_netcdf object.
    R_nc = ROMS_netcdf(args.irradiance_nc_file,Init_ROMS_Irr_Params=(True))

    ## The time step 
    stp = 1

    ## The dataset object from P. Mattern's model out.
    M_nc = Dataset(args.mattern_nc_file) 
    ## Calculating and loading the fields for given wavelength
    irr_field_py = Irradiance_Field_py(R_nc, M_nc, stp)
    ## Making the Eu surface dictionary from irr_field_py output.
    ## This dictionary is the necessary arg. for ocean_color_sol. 
    Eu_surf_py = {}
    for lam in R_nc.wavelengths: 
        Eu_surf_py[lam] = irr_field_py[lam][-1,:,:,2] 
    ocean_color_py = ocean_color_sol(Eu_surf_py, R_nc.Ed0, R_nc.Es0)
    
    ## Visualizing the comparison.
    if args.plot:

        Chlorophyll_Concentration_Comparison(stp, R_nc, M_nc, args.obs_nc_file, ocean_color_py)
     
    
