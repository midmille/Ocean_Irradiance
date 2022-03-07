"""
Created: January 31 2022 
Author: Miles D Miller 

This file follows Shuhba's inquiry from an email on January 21, 2022, 5:56am

"Could we please take a step back, and look at the shapes of the reflectance 
spectra as a function of chlorophyll concentration, for (some of) the phytoplankton types? 
Wavelength on the x-axis, log_10 (reflectance) on the y axis, for a selection of chl values 
(say 0.01, 0.1,1, 10, 100), for each phytoplankton type? if the model is not reproducing 
the right spectral shapes, we cannot expect the inversion to work."
"""

## External Mods
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

## User Mods
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_baird.ocean_irradiance_baird as OIB
import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC
from ocean_irradiance_module import Wavelength_To_RGB


def Solve_Irradiance(PI, N, wavelengths, species, chla_array, z, method='shoot_up'): 
    """
    This function calculates the irradiance field for each species of phytoplankton
    """
    
    ## The Rrs field species dictionary, species names are keys.
    Rrs_field_species ={}

    ## The number of chla profs
    N_chla = chla_array.shape[1]
    ## The number of wavelengths
    N_lam = len(wavelengths)

    ## Looping the species.
    for phy_type in species:
        ## The array in which to store all the irradiance solutions into. 
        Rrs_arr = np.zeros((N_chla, N_lam))
        ## Loop over chla profs 
        for k in range(N_chla):
            ## Looping Wavelengths
            for j,lam in enumerate(wavelengths): 
                
                ## Making a wavelength dependent surface bc
                #if lam == 443: 
                    #PI.Ed0 = 0.1
                    #PI.Es0 = 0.9
                #if lam==551:
                    #PI.Ed0 = 0.9
                    #PI.Es0 = 0.1
                ## The phy object creation
                phy = OI.Phy(z, chla_array[:,k], ESD(phy_type), abscat(lam, phy_type, C2chla='default')[0], abscat(lam, phy_type, C2chla='default')[1])
                ## shoot up method using dutkiewicz model.
                if method == 'shoot_up':
                    ocean_irr_sol = OI.ocean_irradiance_shoot_up(
                                                                 z[0],
                                                                 PI.Ed0,
                                                                 PI.Es0,
                                                                 PI.Euh,
                                                                 abscat(lam, 'water'),
                                                                 PI.coefficients,
                                                                 phy=phy,
                                                                 CDOM=None,
                                                                 N=N,
                                                                 pt1_perc_zbot=True,
                                                                 pt1_perc_phy=True
                                                                 )
                    Eu_surf = ocean_irr_sol[2][-1]
                ## Shubha two stream model.
                if method == 'shubha': 
                    ocean_irr_sol = OIS.ocean_irradiance_shubha(z[0],
                                                                PI.Ed0+PI.Es0, 
                                                                abscat(lam, 'water'), 
                                                                PI.coefficients, 
                                                                phy=phy, 
                                                                CDOM=None, 
                                                                N=N, 
                                                                pt1_perc_zbot =True, 
                                                                pt1_perc_phy=True)
                    Eu_surf = ocean_irr_sol[1][-1]

                Rrs_arr[k, j] = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf) 
     
        ## save irr_field into irr_field_species
        Rrs_field_species[phy_type] = Rrs_arr

    return Rrs_field_species


def Plot_Rrs_Field_Species(Rrs_field_species, wavelengths, species, chlas):
    """
    Plotting the resulting Rrs values
    """

    fig, axes = plt.subplots(nrows=1, ncols=len(species))

    ## Loop species
    for k,phy_type in enumerate(species):
        ax = axes[k]
        Rrs_arr = Rrs_field_species[phy_type]

        for i, chla in enumerate(chlas):
            ax.plot(wavelengths, np.log(Rrs_arr[i, :]), label = chla)

        ax.legend(title='Chla Concentration')
        if k == 0:
            ax.set_ylabel('Log(Rrs)')
        ax.set_xlabel('Wavelength [nm]')
        ax.set_title(phy_type)

    fig.show()

    return 

if __name__ == '__main__':

    ## Params 
    PI = Param_Init()
    
    ## Vertical resolution. 
    N = 200

    ## The values for chla
    chlas = [0.01, 0.1, 1, 10, 100]
    ## Making a full profile for each chla value. 
    chla_array = np.zeros((N, len(chlas)))
    for k, chla in enumerate(chlas):
        chla_array[:,k] = np.full(N,chla)
    ## The accompanying z-grid
    z = np.linspace(-1000, 0, N)

    ## The phytoplankton species
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']

    ## Wavelengths
    wavelengths = PI.wavelengths

    save_file = 'out/Rrs_field_species.p' 
    
    ## Method 
    method = 'shoot_up'

    Rrs_field_species = Solve_Irradiance(PI, N, wavelengths, species, chla_array, z, method=method)
    Plot_Rrs_Field_Species(Rrs_field_species, wavelengths, species, chlas)
