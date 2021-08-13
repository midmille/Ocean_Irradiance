# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 18:01:00 2021

@author: miles

This file implements the ROMS Irradiance resolution study.

Its contents are as follows:
    
    RUN EXPERIMENT [Remote]
    --------------
    
    --> Run this bit from within python code. i.e. create dirbin path. 
    --> Check for output directory of this experiment within dirbin
        -- if directory exist user input on what to do. overwrite? 
    --> Create a array of the independent variable N_irr.
    --> Compile ROMS if told.
    --> Loop N_irr array
        --> Edit N_irr in nemuro.in file. 
            -- replace old instance in file with new.
        --> Edit the name of the roms_his.nc file to reflect N_irr of case. 
            -- This happens in file roms_upwelling.in
            -- replace old instance in file with new.
        --> Run ROMS for these in files.
            -- use subprocess not os for this.
            -- check for non-zero exit code. 
            -- output files
            

NOTE :: Maybe differentiate using flag such that RUN EXPERIMENT is performed only on remote
while the POSTPROCESSING is done locally. This might require saving desired dependent variable
as pickle/np file within working directory. 
            
    RUN POSTPROCESSING [Maybe complete on remote not local]
    -------------------
    
    --> Read in desired output from nc files using the Read_ROMS_Ouutpu module.
    --> Take the highest res to be truth. 
    --> calculate relative difference in given metric.
        -- metric could be: OCx, or Eu_surf. 
    --> If calculations ar costly save output as pickle or np file. 
    
    PLOT [Local]
    ----
    --> Produce visualization of independent N_irr vs dependent metric. 
        
"""
## Global mods
import numpy as np
import matplotlib.pyplot as plt 
import os
import pickle
import sys
import subprocess

## appending ocean _irradiance module to path
ocean_irradiance_module_path = os.path.abspath('../..')
sys.path.append(ocean_irradiance_module_path)

##User mods
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
from ocean_irradiance_module.Ocean_Irradiance_ROMS import Ocean_Irradiance_Field 
from ocean_irradiance_module.Ocean_Irradiance_ROMS import ocean_color_sol
import ocean_irradiance_module.Shell_Script_Tools as SST


def Run_N_irr_Study(N_irr0, roms_hisname0, N_irrs, nc_outdir, run_dirbin, pick_outdir, pick_file_head):
        
        
        ## Checking for existing out directory and making new one
        SST.Make_Out_Dir(nc_outdir)
        SST.Make_Out_Dir(pick_outdir)
        
        
        ## Big Loop over N_irr
        ## --------------------
        ## Setting the old string instance to the current instance.
        old_N_irr_instr = N_irr0
        old_roms_hisname = roms_hisname0
        ## Looop
        for k,N_irr in enumerate(N_irrs): 
            
            ## The new instance will include the current N_irr
            new_N_irr_instr = f'NIRR == {N_irr}'
            new_roms_hisname = f'HISNAME == {nc_outdir}/roms_his_N_irr_{N_irr}.nc'
            
            
            ## First step of the loop is to edit the '.in' files. 
            SST.Edit_ROMS_In_File(f'{run_dirbin}/nemuro.in', old_N_irr_instr, new_N_irr_instr)
            SST.Edit_ROMS_In_File(f'{run_dirbin}/roms_upwelling.in', old_roms_hisname, new_roms_hisname)
            ## RUN ROMS
            ## --------
            out = subprocess.run(['mpirun', '-np', '9', 'romsM', 'roms_upwelling.in'])
            ## Check for completion.
            if out.returncode != 0:
                sys.exit('Non-Zero Retun Code... Exiting ROMS RUN.')
                
            ## netcdf object and pickle save 
            ## Splitting string by white space and grabbing last index.
            nc_file = new_roms_hisname.split()[-1]
            ## Making ROMS_netcdf object
            R_nc = ROMS_netcdf(nc_file,Init_ROMS_Irr_Params=(True))
            ## Saving Object
            pick_path = f'{pick_outdir}/{pick_file_head}{N_irr}.p'
            pickle.dump(R_nc, open(pick_path, 'wb'))
            
            ## Finally setting current instance as old to be replaced in file with next instance.
            old_N_irr_instr = f'NIRR == {N_irr}'
            old_roms_hisname = f'HISNAME == {nc_outdir}/roms_his_N_irr_{N_irr}.nc'
                
        
        ## After the RUN loop change the file back to its originaal state
        ## replace the last used instance in file with the original instance, ie first in list.
        SST.Edit_ROMS_In_File(f'{run_dirbin}/nemuro.in', old_N_irr_instr, N_irr0 )
        SST.Edit_ROMS_In_File(f'{run_dirbin}/roms_upwelling.in', old_roms_hisname, roms_hisname0 )
        
        return 


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs ROMS irradiance N_irr sensitivity study.')
    # parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    parser.add_argument('--run', action='store_true', help="Run Option")
    parser.add_argument('--vis', action='store_true', help="Plot Option")
    args = parser.parse_args()
    
    # file = args.file
    run = args.run 
    plot = args.plot
    
    N_irrs = np.arange(5,40, 2)
    run_dirbin = '~/runs/20210729_test_light'
    out_dir = run_dirbin+'/N_irr_output'
    ## working directory, where the python file is run.
    cwd = os.getcwd()
    ## Where the pickle files will be
    pick_outdir = f'{cwd}/N_irr_pick_out'
    ## The begining of the name of the pickle file. The end will have the N_irr added for differentiation.
    pick_file_head = 'N_irr_'
    
    ## This will be replaced through the replace instance function. Its what is currently written in file.
    N_irr0 = 'NIRR == 16'
    roms_hisname0 = 'HISNAME == output/roms_his.nc'
    
    
    if run:
        print("Running Experiment")
        ## RUN EXPERIMENT
        ## --------------
        Run_N_irr_Study(N_irr0, roms_hisname0, N_irrs, out_dir, run_dirbin, pick_outdir, pick_file_head)
    
    
















































