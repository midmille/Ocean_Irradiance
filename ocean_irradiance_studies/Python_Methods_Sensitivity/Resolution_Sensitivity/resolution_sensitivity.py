"""
Created: October 19 2021 9:26am

Author: Miles Miller
"""

## External Modules
import numpy as np
import matplotlib.pyplot as plt

## User Modules
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR 
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.Wavelength_To_RGB import wavelength_to_rgb
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf


def Run_Study(Ns, R_nc, PI, nstp, N, save_dir, save_files, save_file_head, method):
    """

    """

    ## Number of independent vairables. 
    Nind = len(Ns)
    ## List of irradiance field dictionaries
    irr_fields = []
    for k, N, in enumerate(Ns): 
        
        irr_field = OIR.Irradiance_Run(R_nc, PI, nstp, N, save_dir, save_files[k], method) 
        irr_fields.append(irr_field)

    ## Calculating the truth to be the scipy at the hihest resolution
    irr_field_truth = OIR.Irradiance_Run(R_nc, PI, nstp, N, save_dir, f'{save_file_head}_truth', 'scipy')

    return irr_fields, irr_field_truth


def Process_Results(Ei, irr_fields, irr_field_truth, Ns, wavelengths ):
    """

    """

    ## The number of independent vars.
    Nind = len(Ns)
    ## Calculating the RMS error from the truth.
    abs_rmse = np.zeros(Nind)
    abs_rmse_dict = {}

    for lam in wavelengths: 
        for k, irr_field in enumerate(irr_fields): 
            #abs_rmse[k] = np.sqrt(np.nanmean((irr_field_truth[lam][:,:,:,Ei] - irr_field[lam][:,:,:,Ei])**2)) 
            ## just doing the surface value
            abs_rmse[k] = np.nanmean((irr_field_truth[lam][-1,:,:,Ei] - irr_field[lam][-1,:,:,Ei])) 
        abs_rmse_dict[lam] = abs_rmse

    ## calculating the number of bad points
    num_bad_pts = np.zeros(Nind) 
    num_bad_pts_dict = {}
    for lam in wavelengths: 
        ## Counting bad points
        for k, irr_field in enumerate(irr_fields):
            nlt0 = irr_field[lam][:,:,:,Ei][irr_field[lam][:,:,:,Ei] < 0] 
            ngt1 = irr_field[lam][:,:,:,Ei][irr_field[lam][:,:,:,Ei] > 1]
            num_bad_pts[k] = nlt0.size + ngt1.size
        num_bad_pts_dict[lam] = num_bad_pts


    ## plotting the results
    fig, axes = plt.subplots(1,2)

    x_coords = np.linspace(0, Nind-1, Nind)
    spacing = .1
    for j,lam in enumerate( wavelengths): 
        ## Plotting the mean bias as a function of zbot frac.
        axes[0].plot(Ns, abs_rmse_dict[lam], color=wavelength_to_rgb(lam), label = lam)
        ## Plotting the number of 'bad point' for different zbots
        axes[1].bar(x_coords + j*spacing, num_bad_pts_dict[lam], color=wavelength_to_rgb(lam), label = lam)

        
        
    axes[0].set_title('Resolution Sensitivity Study')
    axes[0].set_ylabel('Mean Bias from Scipy')
    axes[0].set_xlabel('Number of Vertical Layers')
    axes[0].grid()
    axes[0].legend()


    fig.show()

    return abs_rmse_dict, num_bad_pts_dict 
 

    Nind = len(Ns)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Ocean Irradiance ROMS Wrapper')
    parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    parser.add_argument('save_file_head', help = "The name of the pickle file in which the result will be stored" )
    parser.add_argument('save_dir', help = "The name and path of the direcotry in which the result will be stored" )
    parser.add_argument('method', help = "Methods are: Dut, shoot_down, shoot_up, shoot_fp." )
    # parser.add_argument('dest_file', help='Path to Destination Directory. Saved as pickle')
    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    args = parser.parse_args()
    

    ## The independent variable is the number of vertical levels
    Ns = np.arange(10, 300, 15)

    PI = Param_Init()
    R_nc = ROMS_netcdf(args.file)

    ## time step index 
    nstp = 1
 
    Ei = 2

    wavelengths = R_nc.wavelengths

    ## making the file names as functions of N
    save_files = []
    for k, N in enumerate(Ns): 
       save_files.append(f'{args.save_file_head}_{round(N, 3)}.p')
    ## Running the study. 
    irr_fields, irr_field_truth = Run_Study(Ns, R_nc, PI, nstp, N, args.save_dir, save_files, args.save_file_head, args.method)
 
    Process_Results(Ei, irr_fields, irr_field_truth, Ns, wavelengths)
