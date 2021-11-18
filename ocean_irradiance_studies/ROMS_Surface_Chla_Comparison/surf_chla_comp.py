"""
Created: November 10 2021

Author: Miles Miller
"""


## External Mods
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle 
## User made mod
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
import ocean_irradiance_visualization.Plot_Comparison as PC
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat


def Plot_Comparison(chla_ROMS, chla_irrs, species): 
    """
    This plots a comparison between chla from a ROMS output and as calculated by 
    the irradiance model. 
    """ 
    
    fig, ax = plt.subplots()
    
    for k, spec in enumerate(species): 
        ax = PC.Plot_Comparison(ax, chla_ROMS, chla_irrs[:,k], spec, xlim =3, ylim=3)
        #print(species, ' m:', m, ' b:', b)

    ax.set_title('Comparison of Irradiance to ROMS Surface Chla')
    ax.grid()
    ax.set_ylabel('Irradiance Surface Chla')
    ax.set_xlabel('ROMS Surface Chla')
    ax.legend(title='Species')

    fig.show()

    return 

def Run_Diff_Species(species, save_dir, save_files, method, nstp, N):
    """
    This runs the irradiance model over different species.
    """ 

    irr_fields = []
    ## This is the start of the calculations done on the independent variable.                 
    ## Looping over the different bottom depths.
    for k, spec in enumerate(species):
        save_path = f'{save_dir}/{save_files[k]}'
        ## Checking if save file exists
        ## if it doesn't exist then redo calculation
        if os.path.exists(save_path) == False:
            irr_field = {}
            for lam in R_nc.wavelengths:
                print('Current Wavelength:', lam)
                if spec == 'water': 
                    ## if water then  phy abscat must be zero. 
                    phy_abscat= [0,0]
                else:
                    ## Same species for small and large phy. 
                    phy_abscat = abscat(lam, spec)

                irr_field[lam] = OIR.Ocean_Irradiance_Field(
                                                  R_nc.maskr, 
                                                  abscat(lam, 'water'), 
                                                  ## using the same coefficinets for large and small phy.
                                                  phy_abscat, 
                                                  phy_abscat, 
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

        ## A list of the dictionaries for each different zbot.
        irr_fields.append(irr_field) 

    return irr_fields 


if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='Ocean Irradiance ROMS Wrapper')
    parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    parser.add_argument('save_file_name', help = "The name of the pickle file in which the result will be stored" )
    parser.add_argument('save_dir', help = "The name and path of the direcotry in which the result will be stored" )
    parser.add_argument('method', help = "Methods are: Dut, shoot_down, shoot_up, shoot_fp." )
    args = parser.parse_args()

    PI = Param_Init()
    R_nc = ROMS_netcdf(args.file)

    ## settting time step index
    time_step_index = 1
    
    N = 100
    
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn', 'water']

    ## Creating save file name 
    save_files = []
    for spec in species:
        save_files.append(f'{args.save_file_name}_{spec}.p')

    ## Calculating irradiances for different species.
    irr_fields = Run_Diff_Species(species, args.save_dir, save_files, args.method, time_step_index, N)

    ## The reference chla. This is also the input into the irradiance model
    chla_ROMS = R_nc.chl_diatom + R_nc.chl_nanophyt  
    chla_surf_ROMS = np.nanmean(chla_ROMS[time_step_index, -10:, :,:], axis = 0) 
    chla_surf_ROMS = chla_surf_ROMS.flatten()  

    ## storing the chla irr arrays
    chla_irrs = np.zeros((len(chla_surf_ROMS), len(species)))

    for k, spec in enumerate(species): 
 
        ## This is the irradiance model chla. 
        Eu_surf_dict = {}
        for lam in R_nc.wavelengths: 
            Eu_surf_dict[lam] = irr_fields[k][lam][-1, :,:, 2]
    
        chla_irr = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)
        ## flattening the arrays. 
        chla_irrs[:,k] = chla_irr.flatten()

    ## Plotting Comparison. Gonna be a lot of points. 
    Plot_Comparison(chla_surf_ROMS, chla_irrs, species)

    












