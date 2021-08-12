# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 08:47:14 2021

@author: miles

ABOUT 
-----

--> This file seeks to provide a comparison of irradiance profiles with a very strict 
metric. The maximum reltaive difference within the profile between the python code 
and the ROMS FORTRAN implementation. 
--> The worst profiles will then be visualized for a qualitative comparison. 
"""

## Global mods
import numpy as np
import matplotlib.pyplot as plt 
import os
import pickle
import sys

## appending ocean _irradiance module to path
ocean_irradiance_module_path = os.path.abspath('../..')
sys.path.append(ocean_irradiance_module_path)

##User mods
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
from ocean_irradiance_module.Ocean_Irradiance_ROMS import Ocean_Irradiance_Field 
from ocean_irradiance_module.Ocean_Irradiance_ROMS import ocean_color_sol



def Irradiance_Field_py_ROMS(R_nc, nstp):
    """
    

    Parameters
    ----------
    R_nc : ROMS_netcdf object
        Holds all the returned arrays from ROMS output and paramters used in the simulation. 
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
        print('Python irradiance field calculation save file does not exist.')
        print('Redoing calculation... ')
        mask = np.ones((R_nc.nyi, R_nc.nxi))
        irr_field_py = {}
        for lam in R_nc.wavelengths:
            irr_field_py[lam] = Ocean_Irradiance_Field(mask, 
                                              R_nc.ab_wat[lam], 
                                              R_nc.ab_diat[lam], 
                                              R_nc.ab_syn[lam], 
                                              R_nc.chl_diatom[nstp,:,:,:], 
                                              R_nc.chl_nanophyt[nstp,:,:,:], 
                                              R_nc.z_r[nstp,:,:,:], 
                                              R_nc.Ed0, 
                                              R_nc.Es0, 
                                              R_nc.Euh)
        
        pickle.dump(irr_field_py, open(save_path, "wb"))
        print('Python calculation complete and saved')
        
    ## if the save file does exist then just load it. gity 
    elif os.path.exists(save_path) == True:
        print(f'Irradiance save file exists! Loading python calculated irradiance field from file "{save_file}"...')
        Eu_surf_py = pickle.load(open(save_path,'rb'))
        print('Yay, file loaded :)')
        
        
    ## ROMS calculated surface Eu
    ##---------------------------
    irr_field_ROMS = {}
    
    Eu_surf_ROMS[R_nc.wavelengths[0]] = R_nc.Eu1[nstp,:,:,:]
    Eu_surf_ROMS[R_nc.wavelengths[1]] = R_nc.Eu2[nstp,:,:,:]
    
    return Eu_surf_py, Eu_surf_ROMS








if __name__ == '__main__':
    
    file_path = os.getcwd()
    file = f"{file_path}/roms_his_phy.nc"
    R_nc = ROMS_netcdf(file,Init_ROMS_Irr_Params=(True))
    ## The time step 
    nstp = 1
