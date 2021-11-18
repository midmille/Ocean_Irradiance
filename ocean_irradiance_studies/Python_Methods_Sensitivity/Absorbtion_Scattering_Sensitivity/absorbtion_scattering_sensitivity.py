"""
Created: October 21 2021 10:34am 
Author: Miles Miller
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


def Run_Study(a_fracs, b_fracs, R_nc, PI, nstp, N, save_dir, save_files, save_file_head, method): 


    """
    This runs the zbot sensitivity study by calling on the Ocean_Irradiance_ROMS module's 
    Irradiance_Run function.

    Parameters
    ----------
    a_fracs: 1-D Array
        The fraction of the default absorbtion coefficients to be used for irradiance calculation.
    b_fracs: 1-D Array
        The fraction of the default backscattering coefficients to be used for irradiance calculation.
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
                                                      method = method)
            
        ## Saving the truth if calculated.
        pickle.dump(irr_field_truth, open(save_path_truth, 'wb'))
        print('Truth calculation saved')

    ## If the truth pickle file does exist then load it 
    elif os.path.exists(save_path_truth) == True:
        irr_field_truth = pickle.load(open(save_path_truth, 'rb'))


    ## Must be a list because it will contain dictionaris. 
    irr_fields = []
    ## This is the start of the calculations done on the independent variable.                 
    ## Looping over the different bottom depths.
    i=0
    for k, a_frac in enumerate(a_fracs):
        for j, b_frac in enumerate(b_fracs):
            save_path = f'{save_dir}/{save_files[i]}'
            ## Checking if save file exists
            ## if it doesn't exist then redo calculation
            if os.path.exists(save_path) == False:
                irr_field = {}
                for lam in R_nc.wavelengths:
                    print('Current Wavelength:', lam)
                    ab_diat = abscat(lam, 'Diat')
                    ab_diat[0] = ab_diat[0]*a_frac
                    ab_diat[1] = ab_diat[1]*b_frac
                    ab_syn = abscat(lam, 'Syn') 
                    ab_syn[0] = ab_syn[0]*a_frac
                    ab_syn[1] = ab_syn[1]*b_frac

                    irr_field[lam] = OIR.Ocean_Irradiance_Field(
                                                      R_nc.maskr, 
                                                      abscat(lam, 'water'), 
                                                      abscat(lam, 'Diat')*a_frac, 
                                                      abscat(lam, 'Syn') *b_frac, 
                                                      R_nc.chl_diatom[nstp,:,:,:], 
                                                      R_nc.chl_nanophyt[nstp,:,:,:], 
                                                      R_nc.z_r[nstp,:,:,:], 
                                                      PI.Ed0, 
                                                      PI.Es0, 
                                                      PI.Euh,
                                                      PI.coefficients,
                                                      N= N, 
                                                      method = method,
                                                      pt1_perc_zbot = True)  
                           
                ## saving theresult as pickle file
                pickle.dump(irr_field, open(save_path, "wb"))
                print('Python calculation complete and saved')
    
            ## if the save file does exist then just load it. gity 
            elif os.path.exists(save_path) == True:
                #print(f'Irradiance save file exists! Loading python calculated irradiance field from file "{save_files[k]}"...')
                irr_field = pickle.load(open(save_path,'rb'))
                #print('Yay, file loaded :)')
    

            i += 1
            ## A list of the dictionaries for each different zbot.
            irr_fields.append(irr_field) 
      
    return irr_fields, irr_field_truth


def Process_Results(PI, irr_fields, irr_field_truth, a_fracs, b_fracs, wavelengths):
    """

    """

    ## The number absorbtion independent variables.
    N_a = len(a_fracs)
    N_b = len(b_fracs)

    ## Calculating the bias from the truth.
    ## -----------------------------------
    bias_rrs = np.zeros((N_a, N_b))
    bias_rrs_dict = {}
    num_bd_pts = np.zeros((N_a, N_b))
    num_bd_pts_dict = {}

    ## initializng indexing var 
    for lam in wavelengths: 
        i = 0
        for k, a_frac in enumerate(a_fracs): 
            for j, b_frac in enumerate(b_fracs): 
                ## Excluding the bad points from mean. 
                irr_field_truth_bounded = irr_field_truth[lam][-1,:,:,2][~ np.any(np.logical_or(0>irr_fields[i][lam][:,:,:,2], 1<irr_fields[i][lam][:,:,:,2]), axis=0)]
                irr_field_bounded = irr_fields[i][lam][-1,:,:,2][ ~ np.any(np.logical_or(0>irr_fields[i][lam][:,:,:,2], 1<irr_fields[i][lam][:,:,:,2]), axis=0)]
                ## calculating the mean bias 
                bias_rrs[k,j] = np.nanmean((OIR.R_RS(PI.Ed0, PI.Es0, irr_field_truth_bounded) - OIR.R_RS(PI.Ed0, PI.Es0, irr_field_bounded))) 
                ## Calculating the number of bad points.
                num_bd_pts[k,j] = irr_fields[i][lam][-1,:,:,2].size -  irr_field_bounded.size 
                
                i +=1
  
        bias_rrs_dict[lam] = np.copy(bias_rrs)
        num_bd_pts_dict[lam] = np.copy(num_bd_pts)


    ## Calculating the difference in chla
    ## ----------------------------------
    bias_chla = np.zeros((N_a, N_b))
    i=0
    for k, a_frac in enumerate(a_fracs): 
        for j, b_frac in enumerate(b_fracs):
            Eu_surf_dict = {}
            Eu_surf_dict_truth = {}
            for lam in wavelengths:
                ## Excluding the bad points from mean. 
#                irr_field_truth_bounded = irr_field_truth[lam][-1,:,:,2][~ np.any(np.logical_or(0>irr_fields[i][lam][:,:,:,2], 1<irr_fields[i][lam][:,:,:,2]), axis=0)]
#                irr_field_bounded = irr_fields[i][lam][-1,:,:,2][ ~ np.any(np.logical_or(0>irr_fields[i][lam][:,:,:,2], 1<irr_fields[i][lam][:,:,:,2]), axis=0)]

                Eu_surf_dict[lam] = irr_fields[i][lam][-1,:,:,2] 
                Eu_surf_dict_truth[lam] = irr_field_truth[lam][-1,:,:,2]
            bias_chla[k,j] = np.nanmean((OIR.ocean_color_sol(Eu_surf_dict_truth, PI.Ed0, PI.Es0) - OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)))
            i+=1
             

    ## Plotting the results of rrs bias.
    ## --------------------------------
    fig, axes = plt.subplots(len(wavelengths), 2)
    for j,lam in enumerate( wavelengths): 
        Xa_fracs, Yb_fracs = np.meshgrid(a_fracs, b_fracs)
        ## Plotting the mean bias of rrs as a function of a_frac and b_frac.
        im1 = axes[j,0].pcolormesh(Xa_fracs, Yb_fracs, bias_rrs_dict[lam], cmap='nipy_spectral')
        im2 = axes[j,1].pcolormesh(Xa_fracs, Yb_fracs, num_bd_pts_dict[lam], cmap='nipy_spectral')

        fig.colorbar(im1, ax=axes[j,0])
        fig.colorbar(im2, ax=axes[j,1])

        axes[j,0].set_title(f'RRS Mean Bias wavelength = {lam}')
        axes[j,1].set_title(f'Number of Bad Profiles wavelength = {lam}')
        axes[j,0].set_ylabel('Fraction of Scattering')
        axes[j,0].set_xlabel('Fraction of Absorbtion')

    fig.show()


    ## Plotting the results of chla bias.
    ## -----------------------------------
    fig, axes = plt.subplots()
    Xa_fracs, Yb_fracs = np.meshgrid(a_fracs, b_fracs)
    ## Plotting chla bias a a pcolormesh 
    im = axes.pcolormesh(Xa_fracs, Yb_fracs, bias_chla, cmap='nipy_spectral')
    fig.colorbar(im, ax=axes)
    axes.set_title(f'Chla Mean Bias')
    axes.set_ylabel('Fraction of Scattering')
    axes.set_xlabel('Fraction of Absorbtion')

    fig.show()


    return  
 

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
    

    ## The fraction of the default absorbtion coefficients.
    a_fracs = np.linspace(.1, 10, 20)
    ## The fraction of the default scattering coefficients.
    b_fracs = np.linspace(.1, 10, 20)


    PI = Param_Init()
    R_nc = ROMS_netcdf(args.file)

    ## time step index 
    nstp = 1

    ## number of vertical levels in grid
    N = 100

    ## Just the two important wavelengths for now.
    wavelengths = R_nc.wavelengths

    ## making the file names as functions of light_level
    save_files = []
    for k, a_frac in enumerate(a_fracs): 
        for j, b_frac in enumerate(b_fracs):
            save_files.append(f'{args.save_file_head}_{round(a_frac, 3)}_{round(b_frac, 3)}.p')


    ## Running the study. 
    irr_fields, irr_field_truth = Run_Study(a_fracs, b_fracs, R_nc, PI, nstp, N, args.save_dir, save_files, args.save_file_head, args.method)
   
    ## Plotting
    Process_Results(PI, irr_fields, irr_field_truth, a_fracs, b_fracs, wavelengths)


