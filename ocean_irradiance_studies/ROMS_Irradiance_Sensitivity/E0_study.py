# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:56:59 2021

@author: miles

Initial Value Study
-------------------

This study is similiar to the N_irr study in formating but the independent variable to be changed will
be Ed0 and Es0. 

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
    
    ## Two independent variable for this case, Ed0 and Es0, They are dependent of one another. 
    Ed0s = np.arange(.1,.9,.05)
    Es0s = 1 - Ed0s
    
    run_dirbin = '/home/midmille/runs/20210812_wc12'
    ## working directory, where the python file is run.
    cwd = os.getcwd()
    ## Where the pickle files will be
    pick_outdir = f'{cwd}/E0_pick_out'
    ## The begining of the name of the pickle file. The end will have the N_irr added for differentiation.
    pick_file_head = 'E0'
    
    ## This will be replaced through the replace instance function. Its what is currently written in file.
    ## Each independent variable for this experiment must be a string formatted as such
    ##    Ed0  Es0  Euh
    E0_0 = '0.7d0 0.3d0 0.0d0'
    E0s = []
    for Ed0, Es0 in zip(Ed0s, Es0s):
        E0s.append(f'{Ed0}d0 {Es0}d0 0.0d0')
    
    E0_name = 'E0'

    

    if args.run:
        
        print("Running Experiment")
        ## RUN EXPERIMENT
        ## --------------
        RSS.Run(run_dirbin, 'bio_43532.in', 'ocean_43532.in', E0_name, E0s, E0_0, pick_outdir, pick_file_head)
 
    if args.plot: 
        
        ## Use OCx as comparison metric to start.
        max_rel_diff = np.zeros((len(E0s)))
        ## The highest resolution is taken as truth.
        R_nc_true = pickle.load(open(f'{pick_outdir}/{pick_file_head}{E0s[0].split()[0]}.p','rb'))
        OCx_true = R_nc_true.OCx
        ## Looping over the independent variable.
        nstp = 1
        for k,E0 in enumerate(E0s): 
            ## Loading from corresponding pickle file.
            R_nc = pickle.load(open(f'{pick_outdir}/{pick_file_head}{E0.split()[0]}.p','rb'))
            print(E0)
            print('k',k)
            ## Calculating relative difference from truth. 
            max_rel_diff[k] = np.max(abs(OCx_true[nstp,:,:] - R_nc.OCx[nstp,:,:]) / OCx_true[nstp,:,:])
            print(np.max(abs(OCx_true[nstp,:,:] - R_nc.OCx[nstp,:,:]) / OCx_true[nstp,:,:]))
            
        
        ## PLOT
        ## ----
        
        fig,ax = plt.subplots()
        
        ax.plot(Ed0s, max_rel_diff)
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
