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

## appending ocean _irradiance module to path
ocean_irradiance_module_path = os.path.abspath('../..')
sys.path.append(ocean_irradiance_module_path)

##User mods
import Run_Sensitivity_Study as RSS


def main(args):
   
    # file = args.file
    args.run 
    args.plot
    
    N_irrs = np.arange(5,40, 2)
    run_dirbin = '/home/midmille/runs/20210729_test_light'
    ## working directory, where the python file is run.
    cwd = os.getcwd()
    ## Where the pickle files will be
    pick_outdir = f'{cwd}/N_irr_pick_out'
    ## The begining of the name of the pickle file. The end will have the N_irr added for differentiation.
    pick_file_head = 'N_irr_'
    
    ## This will be replaced through the replace instance function. Its what is currently written in file.
    N_irr0 = 16
    

    if args.run:
        
        print("Running Experiment")
        ## RUN EXPERIMENT
        ## --------------
        RSS.Run(run_dirbin, 'nemuro.in', 'NIRR', N_irrs, N_irr0, pick_outdir, pick_file_head)
 
    if args.plot: 
        
        ## Use OCx as comparison metric to start.
        max_rel_diff = np.zeros((len(N_irrs)))
        ## The highest resolution is taken as truth.
        R_nc_true = pickle.load(open(f'{pick_outdir}/{pick_file_head}{N_irrs[-1]}.p','rb'))
        OCx_true = R_nc_true.OCx
        ## Looping over the independent variable.
        nstp = 1
        for k,N_irr in enumerate(N_irrs): 
            ## Loading from corresponding pickle file.
            R_nc = pickle.load(open(f'{pick_outdir}/{pick_file_head}{N_irr}.p','rb'))
            print(N_irr)
            print('k',k)
            ## Calculating relative difference from truth. 
            max_rel_diff[k] = np.max(abs(OCx_true[nstp,:,:] - R_nc.OCx[nstp,:,:]) / OCx_true[nstp,:,:])
            print(np.max(abs(OCx_true[nstp,:,:] - R_nc.OCx[nstp,:,:]) / OCx_true[nstp,:,:]))
            
        
        ## PLOT
        ## ----
        
        fig,ax = plt.subplots()
        
        ax.plot(N_irrs, max_rel_diff)
        ax.grid()
        ax.set_title('N_irr Sensitivity Study')
        ax.set_xlabel('N_irr [Number of Edges in Irradiance Grid]')
        ax.set_ylabel('Relative Error [Highest Resolution == Truth]')
        
        fig.show()
        
        return 


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Runs ROMS irradiance N_irr sensitivity study.')
    # parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    parser.add_argument('--run', action='store_true', help="Run Option")
    parser.add_argument('--plot', action='store_true', help="Plot Option")
    args = parser.parse_args()
    
    main(args)
   
                               















































