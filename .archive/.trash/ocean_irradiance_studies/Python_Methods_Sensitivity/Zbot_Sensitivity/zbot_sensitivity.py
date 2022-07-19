"""
Created October 15 2021 12:21pm
Author Miles Miller
"""


## External Modules
import numpy as np
import os
import pickle 
import matplotlib.pyplot as plt

## User Modules
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR 
import ocean_irradiance_module.Ocean_Irradiance as OI
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.Wavelength_To_RGB import wavelength_to_rgb


def Run_Study(light_levels, R_nc, PI, nstp, N, save_dir, save_files, save_file_head, method): 


    """
    This runs the zbot sensitivity study by calling on the Ocean_Irradiance_ROMS module's 
    Irradiance_Run function.

    Parameters
    ----------

    zbots: 1-D Array
        The array of zbots as a fraction of the zbot at the .1% pure water light level. 
    R_ncs: List of ROMS_netcdf objects
        The z_r within the object will be editied such that the bottom most z coordinate will 
        correspond to the given zbots.
    PI: Param_Init object
        ocean_irradiance_module.PARAMS Param_Init class. 
    nstp: Float
        The index of desired time step of ROMS out nc file.
    save_dir: String
        The path to the out directory for the pickle files of the experiment. 
    save_files: List of Strings
        The name of the pickle file corresponding to the given zbot. 
    save_file_head: String
        The begining of the name of the save files without the corresponding zbot value.
    method: String
        The irradiance solution solver method name to be used.

    """
    ## The file name and path to save the truth into
    save_path_truth = f'{save_dir}/{save_file_head}_truth.p'
    irr_field_truth = {}
    ## truth taken to be scipy method at same resolution and 
    ## zbot taken at default .1% light level. 
    if os.path.exists(save_path_truth) == False: 
        print('calculating the scipy truth')
        for lam in R_nc.wavelengths:
            ## The calculation using scipy
            irr_field_truth[lam] = OIR.Ocean_Irradiance_Field(
                                                      R_nc.maskr, 
                                                      abscat(lam, 'water'), 
                                                      abscat(lam, 'Diat'), 
                                                      abscat(lam, 'Syn'), 
                                                      R_nc.chl_diatom[nstp,:,:,:], 
                                                      R_nc.chl_nanophyt[nstp,:,:,:], 
                                                      R_nc.z_r[nstp,:,:,:], 
                                                      PI.Ed0, 
                                                      PI.Es0, 
                                                      PI.Euh,
                                                      PI.coefficients,
                                                      N= N, 
                                                      method = 'scipy')
            
        ## Saving the truth if calculated.
        pickle.dump(irr_field_truth, open(save_path_truth, 'wb'))
        print('Truth calculation saved')

    ## If the truth pickle file does exist then load it 
    elif os.path.exists(save_path_truth) == True:
        irr_field_truth = pickle.load(open(save_path_truth, 'rb'))


    ## Must be a list because it will contain dictionaris. 
    irr_fields = []
    zbots = np.zeros(len(light_levels)) 
    ## This is the start of the calculations done on the independent variable.                 
    ## Looping over the different bottom depths.
    for k, light_level in enumerate(light_levels):
        save_path = f'{save_dir}/{save_files[k]}'
        ## Checking if save file exists
        ## if it doesn't exist then redo calculation
        if os.path.exists(save_path) == False:
            irr_field = {}
            for lam in R_nc.wavelengths:
                ## calculating the .1% zbot for the given wavelength
                a_wat, b_wat = abscat(lam, 'water')
                c_wat  = (a_wat+b_wat)/(PI.v_d)
                zbot = OI.zbot_func(PI.Ed0, c_wat, light_level)
                ## zbot is a fraction of the default value.
                zbots[k] = zbot
                ## z_r at the designated time index
                zbot_arr = np.copy(R_nc.z_r[nstp,0,:,:])
                ## Where the current zbot is shallower than bathemetry,
                ## zbot becomnes the current given zbot.
                zbot_arr[:,:][zbot > zbot_arr[:,:]] = zbot

                       
                print('Current Wavelength:', lam)
                irr_field[lam] = OIR.Ocean_Irradiance_Field(
                                                  R_nc.maskr, 
                                                  abscat(lam, 'water'), 
                                                  abscat(lam, 'Diat'), 
                                                  abscat(lam, 'Syn'), 
                                                  R_nc.chl_diatom[nstp,:,:,:], 
                                                  R_nc.chl_nanophyt[nstp,:,:,:], 
                                                  R_nc.z_r[nstp,:,:,:], 
                                                  PI.Ed0, 
                                                  PI.Es0, 
                                                  PI.Euh,
                                                  PI.coefficients,
                                                  N= N, 
                                                  method = method,
                                                  zbot_arr = zbot_arr, 
                                                  pt1_perc_zbot = False)  
                       
            ## saving theresult as pickle file
            pickle.dump(irr_field, open(save_path, "wb"))
            print('Python calculation complete and saved')

        ## if the save file does exist then just load it. gity 
        elif os.path.exists(save_path) == True:
            print(f'Irradiance save file exists! Loading python calculated irradiance field from file "{save_files[k]}"...')
            irr_field = pickle.load(open(save_path,'rb'))
            print('Yay, file loaded :)')


        ## A list of the dictionaries for each different zbot.
        irr_fields.append(irr_field) 
      
    return irr_fields, irr_field_truth, zbots


def Process_Results(Ei, irr_fields, irr_field_truth, zbots, light_levels, wavelengths):
    """
    This file calulates the RMS error from the truth for each zbot choice and then plots the result. 

    Parameters
    ----------
    irr_fields: List
        A list of dictionaries correspnding to each different zbot.
    irr_field_truth: Dictionary
        A dictionary of the irradiance fields taken to be the truth. 
    zbots: 1-D Array
        The independent variable of the different zbots used. 

    Returns
    -------

    """

    ## The number of independent vars.
    Nind = len(zbots)
    ## Calculating the RMS error from the truth.
    abs_rmse = np.zeros(Nind)
    abs_rmse_dict = {}

    for lam in wavelengths: 
        for k, irr_field in enumerate(irr_fields): 
            #abs_rmse[k] = np.sqrt(np.nanmean((irr_field_truth[lam][:,:,:,Ei] - irr_field[lam][:,:,:,Ei])**2)) 
            ## just doing the surface value
            abs_rmse[k] = np.nanmean((irr_field_truth[lam][-1,:,:,Ei] - irr_field[lam][-1,:,:,Ei])) 
        abs_rmse_dict[lam] = np.copy(abs_rmse)

    ## calculating the number of bad points
    num_bad_pts = np.zeros(Nind) 
    num_bad_pts_dict = {}
    for lam in wavelengths: 
        ## Counting bad points
        for k, irr_field in enumerate(irr_fields):
            nlt0 = irr_field[lam][:,:,:,Ei][irr_field[lam][:,:,:,Ei] < 0] 
            ngt1 = irr_field[lam][:,:,:,Ei][irr_field[lam][:,:,:,Ei] > 1]
            num_bad_pts[k] = nlt0.size + ngt1.size
        num_bad_pts_dict[lam] = np.copy(num_bad_pts)


    ## plotting the results
    fig, axes = plt.subplots(1,2)

    x_coords = np.linspace(0, Nind-1, Nind)
    spacing = .1
    for j,lam in enumerate( wavelengths): 
        ## Plotting the mean bias as a function of zbot frac.
        axes[0].plot(light_levels, abs_rmse_dict[lam], color=wavelength_to_rgb(lam), label = lam)
        ## Plotting the number of 'bad point' for different zbots
        axes[1].bar(x_coords + j*spacing, num_bad_pts_dict[lam], color=wavelength_to_rgb(lam), label = lam)

        
        
    axes[0].set_title('zbot Sensitivity Study')
    axes[0].set_ylabel('Mean Bias from Scipy')
    axes[0].set_xlabel('Fraction of Ed0')
    axes[0].grid()
    axes[0].legend()


    fig.show()

    return abs_rmse_dict, num_bad_pts_dict 
   

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
    

    ## The fraction of the default .1% zbot.
    light_levels = np.linspace(.001, .5, 20)


    PI = Param_Init()
    R_nc = ROMS_netcdf(args.file)

    ## time step index 
    nstp = 1

    ## number of vertical levels in grid
    N = 100

    ## The irradiance index, only needed for post processing.
    ## Ei = 0 = Ed, Ei = 1 = Es, Ei = 2 = Eu
    Ei = 2

    ## Just the two important wavelengths for now.
    wavelengths = R_nc.wavelengths

    ## making the file names as functions of light_level
    save_files = []
    for k, light_level in enumerate(light_levels): 
       save_files.append(f'{args.save_file_head}_{round(light_level, 3)}.p')


    ## Running the study. 
    irr_fields, irr_field_truth, zbots = Run_Study(light_levels, R_nc, PI, nstp, N, args.save_dir, save_files, args.save_file_head, args.method)
   

    ## Calculating RMSE and plotting. 
    abs_rmse_dict, num_bad_pts_dict = Process_Results(Ei, irr_fields, irr_field_truth, zbots, light_levels, wavelengths)











    

