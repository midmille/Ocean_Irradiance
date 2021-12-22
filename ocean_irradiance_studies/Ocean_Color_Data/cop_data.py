"""
Created: December 21 2021

Author: Miles Miller

This file is to compared the irradiance model to the level 4 chla data from copernicus.
https://resources.marine.copernicus.eu/products

"""


## User Modules 
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC

## External Modules
from netCDF4 import Dataset
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np


def Run_Irr_Comp_Cop(PI, cop_nc, save_dir, save_file_head, wavelengths, N, phy_type, plot=True):
    """
    This file is to calculate the chla using the irradiance model from the chla profiles of the copernicus data.
    Then compare the resulting surface chla to that of copernicus.
    """
    
    ## Copernicus data. 
    ## The first index is time, I only downloaded one step.
    cop_chla = np.flip(cop_nc.variables['chl'][0, :,:,:], axis=0)
    ## Copernicus Z grid is positive downwards.
    cop_z = -np.flip(cop_nc.variables['depth'][:])
    nzi, nyi, nxi = cop_chla.shape
    

    ## The name of the pickle file that the data is being saved into.
    save_path = f'{save_dir}/{save_file_head}_{phy_type}.p' 

    if os.path.exists(save_path):
        ## Load the calculation results.
        irr_field = pickle.load(open(save_path, 'rb'))
        print('Irr field loaded from file')
        
    elif os.path.exists(save_path) == False:
        ## Run the calculation.
        irr_field = {}
        ## Loop wavelengths 
        for lam in wavelengths:
            i=0
            irr_arr = np.zeros((N, nyi, nxi, 4))
            ## Loop the copernicus domain
            for xi in range(nxi):
                for yi in range(nyi):
                    print(f'{i}/{nyi*nxi}')
                    ## Phytoplankton object
                    phy = OI.Phy(cop_z, cop_chla[:,yi,xi], ESD(phy_type), abscat(lam, phy_type)[0], abscat(lam, phy_type)[1])
                    
                    ## Calculating Irradiances. 
                    ocean_irr_sol = OI.ocean_irradiance_shoot_up(
                                                                 cop_z[0], 
                                                                 PI.Ed0, 
                                                                 PI.Es0, 
                                                                 PI.Euh,
                                                                 abscat(lam, 'water'), 
                                                                 PI.coefficients, 
                                                                 phy=phy,
                                                                 CDOM = None, 
                                                                 N=N,
                                                                 pt1_perc_zbot=True,
                                                                 pt1_perc_phy=True,
                                                                 )
                    ## Ed, Es, Eu
                    irr_arr[:, yi, xi, 0] = ocean_irr_sol[0]
                    irr_arr[:, yi, xi, 1] = ocean_irr_sol[1]
                    irr_arr[:, yi, xi, 2] = ocean_irr_sol[2]
                    irr_arr[:, yi, xi, 3] = ocean_irr_sol[3]
                    
                    i+=1

            ## Save irr array into the dict.
            irr_field[lam] = np.copy(irr_arr)
            
        ## Save result into pickle file
        pickle.dump(irr_field, open(save_path, 'wb'))
        print('Irradiance calculation saved')


    ## Now to caluclate chla. 
    ## Making the Eu at surface dict for ocean color function.
    Eu_surf_dict = {}
    for lam in wavelengths:
        Eu_surf = irr_field[lam][-1,:,:,2]
        ## Calculating the Rrs.
        if lam == 443: 
            Rrs_443 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf)
        if lam == 551: 
            Rrs_551 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf)
        Eu_surf_dict[lam] = Eu_surf

    ## Calculating the chla 
    irr_chla = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0) 
 
    ## Plotting the comparison. 
    if plot:
        fig, ax = plt.subplots()
        ax = PC.Plot_Comparison(ax, cop_chla[-1,:,:].flatten(), irr_chla.flatten(), label = None, xlim = 1.5, ylim=1.5)
        ax.set_xlabel('Copernicus Chla')
        ax.set_ylabel('Irradiance Model Chla')
        fig.show()

    return  cop_chla, irr_chla


def Loop_Species_Irr_Comp_Cop(PI, cop_nc, save_dir, save_file_head, wavelengths, N, species): 
    """
    To loop the species for the comparison of the irradiance model to the copernicus data.
    """

    fig, ax = plt.subplots()

    for k, phy_type in enumerate(species):

        cop_chla, irr_chla = Run_Irr_Comp_Cop(PI, cop_nc, save_dir, save_file_head, wavelengths, N, phy_type, plot=False)
        ax = PC.Plot_Comparison(ax, cop_chla[-1,:,:].flatten(), irr_chla.flatten(), label=phy_type, xlim=1.5, ylim=1.5)

    ax.legend()
    fig.show()
  
    return





if __name__ == '__main__':

    import argparse 
    parser = argparse.ArgumentParser(description='Ocean irradiance calculation and comparison to copernicus data.')
    parser.add_argument('nc_file', help='The name of the copernicus nc data file.')
    parser.add_argument('save_dir', help='Directory to save results into')
    parser.add_argument('save_file_head', help='Heading of the file name for the save file.')
    args = parser.parse_args()

    ## Make the file into a Dataset object.
    cop_nc = Dataset(args.nc_file)
 
    ## Params
    PI = Param_Init()
    ## Wavelengths
    wavelengths = [443, 551]
    ## Number of layers for irradiance calculation.
    N = 200
    ## Species of phytoplankton for the a,b coefficients.
    phy_type = 'Diat'
    ## The species that have coefficients.
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']

    ## Running comparison
    #Run_Irr_Comp_Cop(PI, cop_nc, args.save_dir, args.save_file_head, wavelengths, N, phy_type)
    ## Looping over different spcecies
    Loop_Species_Irr_Comp_Cop(PI, cop_nc, args.save_dir, args.save_file_head, wavelengths, N, species)


    
