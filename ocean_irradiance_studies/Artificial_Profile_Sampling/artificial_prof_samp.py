"""
Created: 14 December 2021 
Author: Miles D. Miller

This file is to manage a study that test the sensitivity of the ocean irradiance
algorithm calculated chla to different samplings of the given profile. 
See notes from December 7 2021.
"""


## User made Modules.
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC

## External mods. 
import matplotlib.pyplot as plt
import numpy as np


def Sampling_Sensitivity(PI, zbot, N_irr, phy_type, wavelengths): 
    """

    """

    ## The array of number of sampiling points. 
    N_samps = np.arange(1, 1000, 5) 

    ## Eu Dictionary.
    Eu_surf_dict = {}
    ## Loop wavelengths.
    for lam in wavelengths:
        ## The Eu surface array. 
        Eu_surf = np.zeros(len(N_samps))
        ## Loop over the N_samp as ind. var.
        for k, N_samp in enumerate(N_samps):
            ## Creating the sample grid.
            z_samp = np.linspace(zbot, 0, N_samp)
            phy_samp = OI.artificial_phy_prof(z_samp, 0, 100, .001)
            phy_samp[phy_samp < 1e-6] = 0

            ## Phytoplankton object. 
            phy = OI.Phy(z_samp,
                         phy_samp,
                         ESD(phy_type),
                         abscat(lam, phy_type, C2chla='default')[0],
                         abscat(lam, phy_type, C2chla='default')[1])
    
            ## Calculating the Irradiance
            irr_out = OI.ocean_irradiance_shoot_up(zbot,
                                                   PI.Ed0, 
                                                   PI.Es0, 
                                                   PI.Euh, 
                                                   abscat(lam, 'water'),
                                                   PI.coefficients,
                                                   phy = phy,
                                                   CDOM=None, 
                                                   N = N_irr, 
                                                   pt1_perc_zbot = True, 
                                                   pt1_perc_phy = False)
            ## Eu at surface.
            Eu = irr_out[2]
            Eu_surf[k] = irr_out[2][-1]
        ## The Eu surface dicionary. 
        Eu_surf_dict[lam] = np.copy(Eu_surf)

    ## Calculate chla.
    chla = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)

    ## To calculate the bias as a function of N_samp, the 
    ## truth is taken to be the chla with the largest N_samp. 
    chla_bias = np.zeros(len(N_samps))
    for k, N_samp in  enumerate(N_samps): 
        ## The last chla in the array corresponds to 'truth'. 
        chla_bias[k] = chla[k] - chla[-1] 
    

    return N_samps, chla_bias, Eu_surf_dict, phy_samp, z_samp, Eu, irr_out[3]


def Plot_Sampling_Sensitivity(N_samps, chla_bias): 
    """
    The purpose of this function is to plot the chla bias as a function
    of the number of sampling points used.
    """
 
    fig, ax = plt.subplots()
    
    ax.plot(N_samps, chla_bias)
    ax.grid()
    ax.set_ylabel('The Chla Bias [mg m^-3]')
    ax.set_xlabel('The Number of Sampling Points')
    ax.set_title('Number of Sample Points Sensitivity')

    fig.show()


if __name__ == '__main__': 

    PI = Param_Init()
    zbot = -500
    N_irr = 200
    phy_type = 'Diat'
    wavelengths = [443, 551]

    N_samps, chla_bias, Eu_surf_dict, phy_samp, z_samp, Eu, z_out = Sampling_Sensitivity(PI, zbot, N_irr, phy_type, wavelengths)

    Plot_Sampling_Sensitivity(N_samps, chla_bias)

