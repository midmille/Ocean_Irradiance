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
import ocean_irradiance_module.Wavelength_To_RGB as W2RGB
from ocean_irradiance_module.PARAMS import Param_Init



def Load_Fields(irr_out_dir):

    """
    This function is to load the pickle files of the irradiances resulting from the different methods. 

    """

    
    irr_out_scipy = pickle.load(open(f'ocean_irradiance_out/oi_out_scipy.p', 'rb')) 
    irr_out_shoot_down = pickle.load(open(f'ocean_irradiance_out/oi_out_shootdown.p', 'rb')) 
    irr_out_shoot_up = pickle.load(open(f'ocean_irradiance_out/oi_out_shootup.p', 'rb')) 
    irr_out_dut = pickle.load(open(f'ocean_irradiance_out/oi_out_dutkiewicz.p', 'rb')) 


    return irr_out_scipy, irr_out_shoot_down, irr_out_shoot_up, irr_out_dut



def Plot_Irradiance_Fields(Ei, fields, methods, wavelengths, R_nc, plot_surface=True, plot_min=True, plot_max=True, plot_surf_diff=True ): 

    """ 
    This function is to plot the fields for comparison.
    """


    ## The plotting starts here 
    
    ## The number of rows in the plot
    N_rows = len(wavelengths)
    ## The number of columns in the plot 
    N_cols = len(fields)

    irr_str = ['Ed', 'Es', 'Eu', 'zarr']

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
                    fig = PF.Plot_Fields(N_rows, N_cols, i, field[lam][surf_i,:,:,Ei], R_nc.lat_rho, R_nc.lon_rho, title, f'{irr_str[Ei]} at Surface', vmin=vmin, vmax=vmax)
                else : 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, field[lam][surf_i,:,:,Ei], R_nc.lat_rho, R_nc.lon_rho, title, f'{irr_str[Ei]} at Surface', fig=fig, vmin=vmin, vmax=vmax)

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
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.min(field[lam][:,:,:,Ei], axis=0),  R_nc.lat_rho, R_nc.lon_rho, title, f'Minimum Value of {irr_str[Ei]}', vmin=vmin, vmax=vmax)
                else : 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.min(field[lam][:,:,:,Ei], axis=0), R_nc.lat_rho, R_nc.lon_rho, title, f'Minimum Value of {irr_str[Ei]}', fig=fig, vmin=vmin, vmax=vmax)

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
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.max(field[lam][:,:,:,Ei], axis=0),  R_nc.lat_rho, R_nc.lon_rho, title, f'Maximum Value of {irr_str[Ei]}', vmin=vmin, vmax=vmax)
                else : 
                    fig = PF.Plot_Fields(N_rows, N_cols, i, np.max(field[lam][:,:,:,Ei], axis=0), R_nc.lat_rho, R_nc.lon_rho, title, f'Maximum Value of {irr_str[Ei]}', fig=fig, vmin=vmin, vmax=vmax)

        fig.show()       


    if plot_surf_diff: 
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
                
                ## if j==1 assumes that the scipy solution is the first method 
                ## The scipy solution is taken as truth and thus its magnitude will eb plotted. 
                if j==0:
                    ## The initial figure is created and returned by the Plot_Fields function          
                    if i == 1: 
                        fig = PF.Plot_Fields(N_rows, N_cols, i, field[lam][surf_i,:,:,2],  R_nc.lat_rho, R_nc.lon_rho, title, r'$E_u(z =0)$', vmin=vmin, vmax=vmax)
                    else : 
                        fig = PF.Plot_Fields(N_rows, N_cols, i, field[lam][surf_i,:,:,2], R_nc.lat_rho, R_nc.lon_rho, title, r'$E_u(z =0)$', fig=fig, vmin=vmin, vmax=vmax)
                else : 
                    Eu0_diff = (fields[0][lam][-1,:,:,2] - field[lam][-1,:,:,2]) / fields[0][lam][-1,:,:,2]
                    ## Assumes that the scipy method is the first method.
                    fig = PF.Plot_Fields(N_rows, N_cols, i, Eu0_diff, R_nc.lat_rho, R_nc.lon_rho, title, r'$\frac{Eu0_{\mathrm{scipy}-Eu0{\mathrm{numerical}}}}{Eu0_{\mathrm{scipy}}}$', fig=fig, vmin=vmin, vmax=vmax)




        fig.show()       



    return 


def Count_Bad_Profs(array):

   	
    ## Number of points less than zero 
    mask_lt0 = np.min(array, axis=0) < -1e-7
    nlt0 = array[0,:,:,:][mask_lt0].size
    ## Number of points greater than one 
    mask_gt1 = np.max(array, axis=0) >1
    ngt1 = array[0,:,:,:][mask_gt1].size
        
    return nlt0, ngt1


def Plot_Number_Bad_Profiles(Ei, fields, methods, wavelengths, R_nc):
    """
    This function will count the number of bad profiles in each method at each wave length.
    A bad profile is considered one in which the profile is either negative or greater than one. 
    """

    ## [The total number of points.]
    num_profs = fields[0][443][0, :,:,3][~np.isnan(fields[0][443][0, :,:,3])].size
    ## number of point less than zero for each method.
    Nlt0 = np.zeros((len(fields), len(wavelengths)))
    ## number of points greater than one for each method.
    Ngt1 = np.zeros((len(fields), len(wavelengths)))

    for j, field in enumerate(fields):
        for k, lam in enumerate(wavelengths):
            Nlt0[j,k], Ngt1[j,k] = Count_Bad_Profs(field[lam][:,:,:,:3])
 
    fig, axs = plt.subplots(nrows=1, ncols=len(wavelengths))
    
    pos = [k for k in range(len(methods))]
    for k,lam in enumerate(wavelengths):
        ax = axs[k]
        rgb = W2RGB.wavelength_to_rgb(lam)
        ## number of points greater than one.
        ax.bar(pos, 100*(Ngt1[:,k]/(3*num_profs)), tick_label=methods, color = rgb, align='center', hatch='xx', label=f'{lam} [nm], E > 1 ')
        ## number of points less than 0
        ax.bar(pos, 100*(Nlt0[:,k]/(3*num_profs)), bottom=100*(Ngt1[:,k]/(3*num_profs)), tick_label=methods, color = rgb, align='center', hatch='//', label=f'{lam} [nm], E < 0')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=75)

        if k == 0: 
            ax.set_ylabel('Percentage of ROMS Profiles')
    
        ax.legend() 

    fig.show()


    return 

             
def main(): 
    """
    The main function.

    """

    
    irr_out_dir = '/home/midmille/Ocean_Irradiance/ocean_irradiance_out'
    ROMS_file = '/home/midmille/runs/20210812_wc12/output/wc12_his_43532.nc'
    
    ## The irradiance index [ Edi=0, Esi=1, Eui=2, zarri=2] 
    ## Eu
    Ei = 2
    
    PI = Param_Init()

    #irr_out_scipy, irr_out_shoot_down, irr_out_shoot_up, irr_our_shoot_fp, irr_out_dut = Load_Fields(irr_out_dir)
    fields = Load_Fields(irr_out_dir)
    
    methods = [
               'Scipy',
               'Shoot Down', 
               'Shoot Up', 
               'Dutkiewicz'
              ]

    ## The two wavelengths necessary for OCx chla calculation 
    wavelengths = PI.wavelengths

    ## The ROMS output oobject 
    R_nc = RRO.ROMS_netcdf(ROMS_file)

    Plot_Irradiance_Fields( Ei, fields, methods, wavelengths, R_nc, plot_surface=False, plot_min=False, plot_max=False, plot_surf_diff = True) 
     
    Plot_Number_Bad_Profiles(Ei, fields, methods, wavelengths, R_nc)
   
    return fields

if __name__ == '__main__':

    fields = main() 
