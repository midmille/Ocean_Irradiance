# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:41:40 2021

@author: miles

ABOUT
-----
This file contains the function Run_Sensitivity_Study
This function runs sensitivity experiment in ROMS by editing the desired irradiance 
input files then running the romsM executable for each case. 

"""
## Global Mods
import os
import subprocess
import sys
import pickle

## User Mods
import ocean_irradiance_module.Shell_Script_Tools as SST
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 


def Run(run_dirbin, in_file_name, var_name, var_vals, var_val0, pick_outdir,
        pick_file_head):
         
    ## changing working dir to run_dirbin to run ROMS... maybe dir change should be out of loop.
    cwd = os.getcwd()
    os.chdir(run_dirbin)
    
    ## Checking for existing out directory and making new one
    SST.Make_Out_Dir(pick_outdir)
    
    
    ## Big Loop over N_irr
    ## --------------------
    ## Setting the old string instance to the current instance.
    old_instr = var_val0

    ## Looop
    for k,var_val in enumerate(var_vals): 
        
        ## The new instance will include the current N_irr
        new_instr = f'{var_name} == {var_val}'
        
        
        ## First step of the loop is to edit the '.in' files. 
        SST.Edit_ROMS_In_File(in_file_name, old_instr, new_instr)
        
        ## RUN ROMS
        ## --------

        out = subprocess.run(['mpirun', '-np', '9', 'romsM', 'roms_upwelling.in'])
        ## Check for completion.
        if out.returncode != 0:
            sys.exit('Non-Zero Retun Code... Exiting ROMS RUN.')
        
            
        ## netcdf object and pickle save 
        nc_file = 'output/roms_his.nc'
        ## Making ROMS_netcdf object
        R_nc = ROMS_netcdf(nc_file,Init_ROMS_Irr_Params=(True))
        ## Saving Object
        ## Check for white spaces
        ## if list after split greater than one, concatinate values... 
        ## this will be important for initial value study 
        if type(var_val) == str: 
            var_val_list = var_val.split()
            pick_path = f'{pick_outdir}/{pick_file_head}{var_val_list[0]}.p'
        else:
            pick_path = f'{pick_outdir}/{pick_file_head}{var_val}.p'

        pickle.dump(R_nc, open(pick_path, 'wb'))
        
        ## Finally setting current instance as old to be replaced in file with next instance.
        old_instr = f'{var_name} == {var_val}'
            
    
    ## After the RUN loop change the file back to its originaal state
    ## replace the last used instance in file with the original instance, ie first in list.
    SST.Edit_ROMS_In_File(f'{run_dirbin}/nemuro.in', old_instr, var_val0 )
    
    ## Change back to working directory.
    os.chdir(cwd)
    
    return 

