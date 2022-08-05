"""
Created: October 19 2021 9:26am

Author: Miles Miller
"""

## External Modules
import numpy as np
import matplotlib.pyplot as plt
import os

## User Modules
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR 
from ocean_irradiance_module.PARAMS import Param_Init
import ocean_irradiance_module.Wavelength_To_RGB as W2RGB
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf


def Run_Study(Ns, methods, R_nc, PI, save_dir, save_file_head):
    """

    """
    wavelengths = PI.wavelengths

    irr_method_fields = []
    for j, method in enumerate(methods): 
        irr_fields = []
        ## List of irradiance field dictionaries
        for k, N, in enumerate(Ns): 
            savefile = f'{save_dir}/{save_file_head}_{method}_N{N}.p'
            
            ## [load the file in series]
            if os.path.exists(savefile) == True: 
                irr_field = OIR.Irradiance_Run(R_nc, PI, 1, N, savefile, method, wavelengths)
                irr_fields.append(irr_field)

            ## [run in parallel]
            elif os.path.exists(savefile) == False: 
                irr_field = OIR.Irradiance_Run_Parallel(R_nc, PI, 1, N, savefile, method, wavelengths)
                irr_fields.append(irr_field)

        irr_method_fields.append(irr_fields)

    return irr_method_fields


def Count_Bad_Profs(array):

   	
    ## Number of points less than zero 
    mask_lt0 = np.min(array, axis=0) < -1e-7
    nlt0 = array[0,:,:,:][mask_lt0].size
    ## Number of points greater than one 
    mask_gt1 = np.max(array, axis=0) >1
    ngt1 = array[0,:,:,:][mask_gt1].size
        
    return nlt0, ngt1


def Plot_Bad_Profiles(Ns, methods, wavelengths, irr_method_fields): 
    """
    """
    fig, axes = plt.subplots(nrows=len(methods), ncols = 1)

    for j, method in enumerate(methods): 

            
        irr_fields = irr_method_fields[j]
        ax = axes[j]
            
        ## [The total number of points.]
        num_profs = irr_fields[0][443][0, :,:,3][~np.isnan(irr_fields[0][443][0, :,:,3])].size

        Nlt0 = np.zeros((len(irr_fields), len(wavelengths)))
        Ngt1 = np.zeros((len(irr_fields), len(wavelengths)))

        for i, field in enumerate(irr_fields): 
            for k, lam in enumerate(wavelengths): 
                Nlt0[j,k], Ngt1[j,k] = Count_Bad_Profs(field[lam][:,:,:,:3])

        width = 1
        pos = Ns - (len(wavelengths)*width)/2
        for k, lam in enumerate(wavelengths): 
            rgb = W2RGB.wavelength_to_rgb(lam)
            ax.bar(pos, 100*(Ngt1[:,k]/(3*num_profs)), color = rgb, align='center', hatch='xx', label=f'{lam} [nm], E > 1 ', width =width)
            
            ax.bar(pos, 100*(Nlt0[:,k]/(3*num_profs)), bottom=100*(Ngt1[:,k]/(3*num_profs)), color = rgb, align='center', hatch='//', label=f'{lam} [nm], E < 0', width = width)

            pos = pos+width

        ax.set_title(f"{method}")
        

    fig.show()
            

    return


if __name__ == '__main__': 

    PI = Param_Init()

    save_dir ='irr_out' 

    save_file_head = 'irr'

    roms_file = '/home/midmille/runs/20210812_wc12/output/wc12_his_43532.nc'

    wavelengths = PI.wavelengths

    R_nc = ROMS_netcdf(roms_file)

    Ns = np.arange(10, 110, 20)

    methods = ['shootup', 'shootdown', 'scipy', 'dutkiewicz']

    irr_method_fields = Run_Study(Ns, methods, R_nc, PI, save_dir, save_file_head)

    Plot_Bad_Profiles(Ns, methods, wavelengths, irr_method_fields)
