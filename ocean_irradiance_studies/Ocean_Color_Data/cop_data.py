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
from scipy import interpolate


def Run_Irr_Comp_Cop(PI, cop_chla_nc, cop_rrs443_nc, cop_rrs560_nc, save_dir, save_file_head, wavelengths, N, phy_type, plot=True):
    """
    This file is to calculate the chla using the irradiance model from the chla profiles of the copernicus data.
    Then compare the resulting surface chla to that of copernicus.
    """
    
    ## Copernicus data. 
    ## The first index is time, I only downloaded one step.
    cop_chla = np.flip(cop_chla_nc.variables['chl'][0, :,:,:], axis=0)
    ## Copernicus Z grid is positive downwards.
    cop_z = -np.flip(cop_chla_nc.variables['depth'][:])
    nzi, nyi, nxi = cop_chla.shape

    ## Interpolating the copernicus rrs to be on copernicus chla grid.
    ## Chla grid 
    cop_chla_lat = cop_chla_nc.variables['latitude'][:]
    cop_chla_lon = cop_chla_nc.variables['longitude'][:]
    ## rrs grid, only gonna use 443, 551 is same grid.
    cop_rrs_lat = cop_rrs443_nc.variables['lat'][:]
    cop_rrs_lon = cop_rrs443_nc.variables['lon'][:]
    ## the values of rrs.
    cop_rrs443 = cop_rrs443_nc.variables['RRS443'][:]
    cop_rrs560 = cop_rrs560_nc.variables['RRS560'][:]
    ## The actual 2-D interpolation results in a function.
    f_interp_rrs443 = interpolate.interp2d(cop_rrs_lat, cop_rrs_lon, cop_rrs443)
    f_interp_rrs560 = interpolate.interp2d(cop_rrs_lat, cop_rrs_lon, cop_rrs560)
    ## The interpolated values onto the chla grid.
    cop_rrs443 = f_interp_rrs443(cop_chla_lat, cop_chla_lon)
    cop_rrs560 = f_interp_rrs560(cop_chla_lat, cop_chla_lon)

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
            irr_rrs443 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf)
        if lam == 551: 
            irr_rrs551 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf)
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

    return  cop_chla, cop_rrs443, cop_rrs560, irr_chla, irr_rrs443, irr_rrs551


def Loop_Species_Irr_Comp_Cop(PI, cop_chla_nc, cop_rrs443_nc, cop_rrs560_nc, save_dir, save_file_head, wavelengths, N, species): 
    """
    To loop the species for the comparison of the irradiance model to the copernicus data.
    """

    fig, axes = plt.subplots(nrows=1, ncols=3)
    ax1, ax2, ax3 = axes
    for k, phy_type in enumerate(species):

        cop_chla, cop_rrs443, cop_rrs560, irr_chla, irr_rrs443, irr_rrs551 = Run_Irr_Comp_Cop(PI, cop_chla_nc, cop_rrs443_nc, cop_rrs560_nc, save_dir, save_file_head, wavelengths, N, phy_type, plot=False)
        ax1 = PC.Plot_Comparison(ax1, cop_chla[-1,:,:].flatten(), irr_chla.flatten(), label=phy_type, xlim=1.5, ylim=1.5)
        ax2 = PC.Plot_Comparison(ax2, cop_rrs443.flatten(), irr_rrs443.flatten(), label=phy_type, xlim=.02, ylim=.02)
        ax3 = PC.Plot_Comparison(ax3, cop_rrs560.flatten(), irr_rrs551.flatten(), label=phy_type, xlim=.02, ylim=.02)

    ax1.set_ylabel('Irr Model Chla [mg chl-a m^-3]')
    ax2.set_ylabel('Irr Model Rrs [sr^-1]')
    ax3.set_ylabel('Irr Model Rrs [sr^-1]')
    ax1.set_xlabel('Copernicus Chla [mg chl-a m^-3]')
    ax2.set_xlabel('Copernicus Rrs [sr^-1]') 
    ax3.set_xlabel('Copernicus Rrs [sr^-1]') 
    ax1.set_title('Chla')
    ax2.set_title('rrs443')
    ax3.set_title('rrs551')
    ax1.legend()
    ax2.legend()
    ax3.legend()
    fig.show()
  
    return





if __name__ == '__main__':

    import argparse 
    parser = argparse.ArgumentParser(description='Ocean irradiance calculation and comparison to copernicus data.')
    parser.add_argument('chla_nc_file', help='The name of the copernicus nc data file.')
    parser.add_argument('rrs443_nc_file', help='The name of the copernicus nc data file.')
    parser.add_argument('rrs560_nc_file', help='The name of the copernicus nc data file.')
    parser.add_argument('save_dir', help='Directory to save results into')
    parser.add_argument('save_file_head', help='Heading of the file name for the save file.')
    args = parser.parse_args()

    ## Make the file into a Dataset object.
    cop_chla_nc = Dataset(args.chla_nc_file)
    cop_rrs443_nc = Dataset(args.rrs443_nc_file)
    cop_rrs560_nc = Dataset(args.rrs560_nc_file)
 
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
    #Run_Irr_Comp_Cop(PI, cop_chla_nc, args.save_dir, args.save_file_head, wavelengths, N, phy_type)
    ## Looping over different spcecies
    Loop_Species_Irr_Comp_Cop(PI, cop_chla_nc, cop_rrs443_nc, cop_rrs560_nc, args.save_dir, args.save_file_head, wavelengths, N, species)


    
