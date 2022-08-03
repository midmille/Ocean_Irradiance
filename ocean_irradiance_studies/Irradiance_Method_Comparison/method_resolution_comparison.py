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
from ocean_irradiance_module.Wavelength_To_RGB import wavelength_to_rgb
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


if __name__ == '__main__': 

    PI = Param_Init()

    save_dir ='irr_out' 

    save_file_head = 'irr'

    roms_file = '/home/midmille/runs/20210812_wc12/output/wc12_his_43532.nc'

    R_nc = ROMS_netcdf(roms_file)

    Ns = np.arange(10, 200, 20)

    methods = ['shootup', 'shootdown', 'scipy', 'dutkiewicz']

    Run_Study(Ns, methods, R_nc, PI, save_dir, save_file_head)


