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
                                              R_nc.chl_diatom[nstp,:,:,:], 
                                              R_nc.chl_nanophyt[nstp,:,:,:], 
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
        
        
    ## ROMS calculated surface Eu
    ##---------------------------
    irr_field_ROMS = {}

    ## The first wavelength, will need to edit so as to loop over them. 
    irr_field_arr = np.zeros((R_nc.N_irr, R_nc.nyi, R_nc.nxi,4))      
    irr_field_arr[:,:,:,0] = R_nc.Ed1[nstp,:,:,:]
    irr_field_arr[:,:,:,1] = R_nc.Es1[nstp,:,:,:] 
    irr_field_arr[:,:,:,2] = R_nc.Eu1[nstp,:,:,:]
    irr_field_arr[:,:,:,3] = R_nc.z_irr1[nstp,:,:,:]
    irr_field_ROMS[R_nc.wavelengths[0]] = irr_field_arr 

    ## The second wavelength. 
    irr_field_arr = np.zeros((R_nc.N_irr, R_nc.nyi, R_nc.nxi,4))      
    irr_field_arr[:,:,:,0] = R_nc.Ed2[nstp,:,:,:]
    irr_field_arr[:,:,:,1] = R_nc.Es2[nstp,:,:,:] 
    irr_field_arr[:,:,:,2] = R_nc.Eu2[nstp,:,:,:]
    irr_field_arr[:,:,:,3] = R_nc.z_irr2[nstp,:,:,:]
    irr_field_ROMS[R_nc.wavelengths[1]] = irr_field_arr 



    return irr_field_py, irr_field_ROMS


def Irradiance_Field_Diff(nstp, R_nc, irr_index, lam, irr_field_py, irr_field_ROMS, 
                              abs_or_rel='abs', plot=False):
        
        """
        

        Parameters
        ----------
        nstp : Int
            Time Step index. 
        R_nc : ROMS_netcdf object.
            All paramteres and necessary ROMs outputs. 
        irr_index : Int
            The irradiance field to calculate the difference for.
        lam : Integer
            The wavelength to calculate this for, it will serve as the key to the dictionary args.
        irr_field_py : Dictionary
            Keys are the wavelengths, the values are the irradiance arrays indexed as follows:
                (z,y,x,irradiance_field_index). 
        irr_field_ROMS : Dictionary
            Keys are the wavelengths, the values are the irradiance arrays indexed as follows:
                (z,y,x,irradiance_field_index). 
        abs_or_rel: String, OPTIONS=['abs', 'rel']
            Absolute or Relative Error Calculation. 
        plot : Boolean, optional
            True or False. The default is False.

        Returns
        -------
        irr_field_max_diff : TYPE
            DESCRIPTION.

        """
        ## Edititng so that the border of the python code is the same as ROMS
        for lam in R_nc.wavelengths:
            irr_field_py[lam][:,0,:,:] = 0  
            irr_field_py[lam][:,-1,:,:] = 0 
            irr_field_py[lam][:,:,0,:] = 0 
            irr_field_py[lam][:,:,-1,:] = 0 
        
        ## This is to be made into a function that finds the worst difference in profiles

        ## difference, rel or abs depending on flag. 
        if abs_or_rel == 'rel': 
            irr_field_diff = abs((irr_field_py[lam] - irr_field_ROMS[lam])/irr_field_py[lam])
        elif abs_or_rel == 'abs':
            irr_field_diff = abs((irr_field_py[lam] - irr_field_ROMS[lam]))
        else:
            sys.exit('MAYDAY, incorrect flag')
    
        ## The maximum error in each profile, could be relative or absolute. 
        irr_field_max_diff = np.zeros((R_nc.nyi, R_nc.nxi))
          
        
        ## looping over horixontal domain and finding maximum differences in each profile
        for j in range(R_nc.nyi): 
            for i in range(R_nc.nxi):
                irr_field_max_diff[j,i] = np.max(irr_field_diff[:, j, i, irr_index])
        
        ## Taking the mask into accoutn to plot easier
        irr_field_max_diff[R_nc.maskr==0] = np.NaN
        if plot: 
            
            fig,ax = plt.subplots()
            im = ax.pcolormesh(irr_field_max_diff)
            fig.colorbar(im, ax = ax, label=r'max_z($\frac{\mathrm{irrpy - EuROMS}}{Eupy}$)')
            if abs_or_rel == 'rel':
                ax.set_title('Relative Difference in Python and Fortran \n Implementation of Irradiance Code \n Using Maximum Vertical Profile Deviation')
            elif abs_or_rel == 'abs':
                ax.set_title('Absolute Difference in Python and Fortran \n Implementation of Irradiance Code \n Using Maximum Vertical Profile Deviation')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            
            fig.show()
        
        return irr_field_max_diff


if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Ocean Irradiance ROMS Wrapper')
    parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    # parser.add_argument('dest_file', help='Path to Destination Directory. Saved as pickle')
    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    args = parser.parse_args()
    
    file = args.file

    R_nc = ROMS_netcdf(file,Init_ROMS_Irr_Params=(True))
    ## The time step 
    nstp = 1
    ## This is the Es index. 
    irr_index = 2
    ## The wavelength choice for this. 
    lam = R_nc.wavelengths[0]
    
    ## Calculating and loading the fields
    irr_field_py, irr_field_ROMS = Irradiance_Field_py_ROMS(R_nc, nstp)
    
    ## Finding maximum vertical profile difference and plotting as color map. 
    Irradiance_Field_Diff(nstp, R_nc, irr_index, lam, irr_field_py, irr_field_ROMS, 
                              abs_or_rel='abs', plot=True)
    
    
    
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
