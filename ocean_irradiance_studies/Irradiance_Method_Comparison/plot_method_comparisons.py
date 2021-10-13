"""
Created on Wed October 13 11:59:28 2021

@author: Miles Miller
"""

## External Mods
import matplotlib.pyplot as plt
import pickle
import numpy as np
## User Mods
import ocean_irradiance_module.Read_ROMS_Out as RRO
import ocean_irradiance_visualization.Plot_Field as PF 



def Load_Fields(irr_out_dir):

    """
    This function is to load the pickle files of the irradiances resulting from the different methods. 

    """

    
    irr_out_scipy = pickle.load(open(f'{irr_out_dir}/irradiance_out_scipy.p', 'rb')) 
    irr_out_shoot_down = pickle.load(open(f'{irr_out_dir}/irradiance_out_shoot_down_log.p', 'rb')) 
    irr_out_shoot_up = pickle.load(open(f'{irr_out_dir}/irradiance_out_shoot_up_log.p', 'rb')) 
    irr_out_shoot_fp = pickle.load(open(f'{irr_out_dir}/irradiance_out_shoot_fp_20_shots.p', 'rb')) 
    irr_out_dut = pickle.load(open(f'{irr_out_dir}/irradiance_out_dut.p', 'rb')) 


    return irr_out_scipy, irr_out_shoot_down, irr_out_shoot_up, irr_out_shoot_fp, irr_out_dut



def Plot_Irradiance_Fields( plot_surface=True, plot_min=True, plot_max=True): 

    """ 
    This function is to plot the fields for comparison.
    """

    irr_out_dir = '/home/midmille/Ocean_Irradiance/ocean_irradiance_out'
    ROMS_file = '/home/midmille/runs/20210812_wc12/output/wc12_his_43532.nc'
    
    ## The irradiance index [ Edi=0, Esi=1, Eui=2, zarri=2] 
    ## Eu
    Ei = 2 
    

    #irr_out_scipy, irr_out_shoot_down, irr_out_shoot_up, irr_our_shoot_fp, irr_out_dut = Load_Fields(irr_out_dir)
    fields = Load_Fields(irr_out_dir)
    
    methods = [
               'Scipy',
               'Shoot Down', 
               'Shoot Up', 
               'Shoot Fit Point', 
               'Dutkiewicz et al. (2015)'
              ]

    ## The two wavelengths necessary for OCx chla calculation 
    wavelengths = [443, 551] 

    ## The ROMS output oobject 
    R_nc = RRO.ROMS_netcdf(ROMS_file)


    ## The plotting starts here 
    
    ## The number of rows in the plot
    N_rows = len(wavelengths)
    ## The number of columns in the plot 
    N_cols = len(fields)

    ## Looping over the wave lengths and the irradiance fields
    ## The index counter for subplot_loc, index starts at 1
    if plot_surface: 
        vmin = 0
        vmax = .19
        i = 0
        for k, lam in enumerate(wavelengths): 
            for j, field in enumerate(fields): 
                ## The index starts at one in the upper left corner of figure
                i += 1    
            
                ## The title fo each axxes instance in the figure.  
                title = f'{methods[j]} \n' + r' [$\lambda = $' + f'{lam}]'
            

                ## The surface index is -1 for all except the Dut. method for which it is 0.
                ## Assumes that the Dut. field is the last one.
                if j == N_cols-1:
                    surf_i = 0
                else: 
                    surf_i = -1

                ## The initial figure is created and returned by the Plot_Fields function          
                if i == 1: 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, field[lam][surf_i,:,:,Ei], R_nc.lat_rho, R_nc.lon_rho, title, 'Eu at Surface', vmin=vmin, vmax=vmax)
                else : 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, field[lam][surf_i,:,:,Ei], R_nc.lat_rho, R_nc.lon_rho, title, 'Eu at Surface', fig=fig, vmin=vmin, vmax=vmax)

        fig.show()


    if plot_min: 
        vmin = None
        vmax = None
        i = 0
        for k, lam in enumerate(wavelengths): 
            for j, field in enumerate(fields): 
                ## The index starts at one in the upper left corner of figure
                i += 1    
            
                ## The title fo each axxes instance in the figure.  
                title = f'{methods[j]} \n' + r' [$\lambda = $' + f'{lam}]'
            

                ## The surface index is -1 for all except the Dut. method for which it is 0.
                ## Assumes that the Dut. field is the last one.
                if j == N_cols-1:
                    surf_i = 0
                else: 
                    surf_i = -1

                ## The initial figure is created and returned by the Plot_Fields function          
                if i == 1: 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.min(field[lam][:,:,:,Ei], axis=0),  R_nc.lat_rho, R_nc.lon_rho, title, 'Minimum Value of Eu', vmin=vmin, vmax=vmax)
                else : 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.min(field[lam][:,:,:,Ei], axis=0), R_nc.lat_rho, R_nc.lon_rho, title, 'Minimum Value of Eu', fig=fig, vmin=vmin, vmax=vmax)

        fig.show()       


    if plot_max: 
        vmin = None
        vmax = None
        i = 0
        for k, lam in enumerate(wavelengths): 
            for j, field in enumerate(fields): 
                ## The index starts at one in the upper left corner of figure
                i += 1    
            
                ## The title fo each axxes instance in the figure.  
                title = f'{methods[j]} \n' + r' [$\lambda = $' + f'{lam}]'
            

                ## The surface index is -1 for all except the Dut. method for which it is 0.
                ## Assumes that the Dut. field is the last one.
                if j == N_cols-1:
                    surf_i = 0
                else: 
                    surf_i = -1

                ## The initial figure is created and returned by the Plot_Fields function          
                if i == 1: 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.max(field[lam][:,:,:,Ei], axis=0),  R_nc.lat_rho, R_nc.lon_rho, title, 'Maximum Value of Eu', vmin=vmin, vmax=vmax)
                else : 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.max(field[lam][:,:,:,Ei], axis=0), R_nc.lat_rho, R_nc.lon_rho, title, 'Maximum Value of Eu', fig=fig, vmin=vmin, vmax=vmax)

        fig.show()       


    return 
             


