# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 14:59:16 2021

@author: miles

ABOUT
-----

This script is to compare the ROMS calculation of Eu at the surfce and compare that 
to the python calculated Eu at surface. 
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


def Eu_Surface_py_ROMS(R_nc, nstp):
    """
    

    Parameters
    ----------
    R_nc : ROMS_netcdf object
        Holds all the returned arrays from ROMS output and paramters used in the simulation. 
    nstp : Float
        The given time step to calculate for.

    Returns
    -------
    Eu_surf_py : Dictionary
        The keys are the wavelengths and the values are the Eu at surface arrays.
        This is calculated by python code.
    Eu_surf_ROMS : Dictionary
        The keys are the wavelengths and the values are the Eu at surface arrays. 
        This is calculated from within the ROMS Fortran code.
    """
    
    ## The name of the file that the python Eu dict will be saved to as pickle.
    save_file = f'Eu_surf_dict_nstp_{nstp}.p'
    save_dir = 'Eu_surf_out'
    save_path = f'{save_dir}/{save_file}'
    ## Python calculated Eu at surface 
    ##---------------------------------
    ## Checking if save file exists
    ## if it doesn't exist then redo calculation
    if os.path.exists(save_path) == False:
        print('Python Eu at surface calculation save file does not exist.')
        print('Redoing calculation... ')
        mask = np.ones((R_nc.nyi, R_nc.nxi))
        Eu_surf_py = {}
        for lam in R_nc.wavelengths:
            print('Wavelength', lam)
            Eu_surf_py[lam] = Ocean_Irradiance_Field(mask, 
                                              R_nc.ab_wat[lam], 
                                              R_nc.ab_diat[lam], 
                                              R_nc.ab_syn[lam], 
                                              R_nc.chl_diatom[nstp,:,:,:], 
                                              R_nc.chl_nanophyt[nstp,:,:,:], 
                                              R_nc.z_r[nstp,:,:,:], 
                                              R_nc.Ed0, 
                                              R_nc.Es0, 
                                              R_nc.Euh,
                                              N= R_nc.N_irr)[-1,:,:,2]
        
        pickle.dump(Eu_surf_py, open(save_path, "wb"))
        print('Python calculation complete and saved')
        
    ## if the save file does exist then just load it. gity 
    elif os.path.exists(save_path) == True:
        print(f'Eu save file exists! Loading python calculated Eu from file "{save_file}"...')
        Eu_surf_py = pickle.load(open(save_path,'rb'))
        print('Yay, file loaded :)')
        
        
    ## ROMS calculated surface Eu
    ##---------------------------
    Eu_surf_ROMS = {} 
    Eu_surf_ROMS[R_nc.wavelengths[0]] = R_nc.Eu1[nstp,-1,:,:]
    Eu_surf_ROMS[R_nc.wavelengths[1]] = R_nc.Eu2[nstp,-1,:,:]
    
    return Eu_surf_py, Eu_surf_ROMS


def Eu_Rel_Difference(lam, Eu_surf_py, Eu_surf_ROMS):
    """
    Calculates the relative difference between Eu at surface of python code and Fortran 
    ROMS code. 

    Parameters
    ----------
    lam : Int
        The wavelength to calculate this for. 
    Eu_surf_py : Dictionary
        The keys are the wavelengths and the values are the Eu at surface arrays.
        This is calculated by python code.
    Eu_surf_ROMS : Dictionary
        The keys are the wavelengths and the values are the Eu at surface arrays. 
        This is calculated from within the ROMS Fortran code.

    Returns
    -------
    rel_diff : 2-D Array
        The relative difference of Eu at surface of python and fortran code. 

    """
    
    diff = Eu_surf_py[lam] - Eu_surf_ROMS[lam]
    rel_diff = diff / Eu_surf_py[lam]

    
    return diff

    
def Eu_Surf_py_ROMS_comparison(nstp, R_nc, plot=False): 
        
        """
        Main file that performs comparison of Eu at surface for python code an dFortran code.

        Parameters
        ----------
        nstp : Int
            Time step index for calculation.
        R_nc : ROMS_netcdf Object.
            Contains parameters and ROMS output. 

        Returns
        -------
        Eu_surf_py : Dictionary
            The keys for this dictionary are the wavelengths and the values are the 
            Eu surface value arrays for the python code. 
        Eu_surf_ROMS: Dictionary
            The keys for this dictionary are the wavelengths and the values are the 
            Eu surface value arrays for the python code. 

        """
        
        ## Calculating respective Eu at surfaces for both code frameworks
        ## saves or loads from file.
        
        Eu_surf_py, Eu_surf_ROMS = Eu_Surface_py_ROMS(R_nc, nstp)
        
        ## Chosen wavelength
        lam = R_nc.wavelengths[1]
        ## Editing Eu_surf_py such that its boundary is also 0 like ROMS
        for lam in R_nc.wavelengths:
            Eu_surf_py[lam][0,:] = 0  
            Eu_surf_py[lam][-1,:] = 0  
            Eu_surf_py[lam][:,0] = 0  
            Eu_surf_py[lam][:,-1] = 0         
        ## Accounting for mask
        Eu_surf_py[lam][R_nc.maskr == 0] = np.NaN
        Eu_surf_ROMS[lam][R_nc.maskr == 0] = np.NaN
            
        ## Finding the relative error for given wavelength 
        rel_diff = Eu_Rel_Difference(lam, Eu_surf_py, Eu_surf_ROMS)
        
        ##  Plotting the relative difference
        if plot:
            fig,ax = plt.subplots()
            im = ax.pcolormesh(rel_diff)
            fig.colorbar(im, ax = ax, label=r'$\frac{\mathrm{Eupy - EuROMS}}{Eupy}$')
            ax.set_title('Relative Difference in Python and Fortran \n Implementation of Irradiance Code \n Using Eu at Surface Metric')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
           
            fig.show()
        
        
        return Eu_surf_py, Eu_surf_ROMS 


def Compare_OCx(nstp, R_nc, plot=False):
    
    ## Ocean_Color 
    ##--------------
    ocean_color_py = ocean_color_sol(Eu_surf_py, R_nc.Ed0, R_nc.Es0)
    # ocean_color_ROMS = ocean_color_sol(Eu_surf_ROMS, R_nc.Ed0, R_nc.Es0)
    ocean_color_ROMS = R_nc.OCx[nstp,:,:]
    ocean_color_ROMS[R_nc.maskr == 0] = np.NaN
    
    if plot: 

        fig,(ax1,ax2) = plt.subplots(1,2)
        #vmax = max(np.max(ocean_color_py), np.max(ocean_color_ROMS))
        #vmin = min(np.min(ocean_color_py), np.min(ocean_color_ROMS))
        vmin,vmax = (0,.5)
        cmap = 'nipy_spectral'
             
        im1 = ax1.pcolormesh(ocean_color_py,vmin=vmin, vmax=vmax, cmap=cmap)
        im2 = ax2.pcolormesh(ocean_color_ROMS,vmin=vmin, vmax=vmax, cmap=cmap)
        
        fig.colorbar(im1, ax = [ax1,ax2] , label=r'chl_a')
        #fig.colorbar(im2, ax = ax2, label=r'chl_a')
    
        ax1.set_title('Ocean Color Python')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax2.set_title('Ocean Color ROMS')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
 
        fig.show()
    
    return ocean_color_py, ocean_color_ROMS

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
    nstp = 3
    
    ## Ucomment this little section to run over all time steps, each step will be saved
    ## as a seperate pickle file. 
    # nstps = np.shape(R_nc.Ed1)[0]
    ##Looping over the time steps
    # for nstp in range(nstps):
    #     Eu_surf_py, Eu_surf_ROMS = Eu_Surface_py_ROMS(R_nc, nstp)
    
    Eu_surf_py, Eu_surf_ROMS = Eu_Surf_py_ROMS_comparison(nstp, R_nc, plot=args.plot)
    
    ocean_color_py, ocean_color_ROMS = Compare_OCx(nstp, R_nc, plot=args.plot)

    
    
    

