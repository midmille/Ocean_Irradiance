"""
Created: January 31 2021
Author: Miles D Miller 


This file follows Shuhba's inquiry from an email on January 21, 2022, 5:56am

"It would also be good to plot on a log-log scale chl on the x-axis and 
Rrs ratio (say Rrs443/Rrs555) or the closest wavelength pair that you have,
again for individual phytoplankton types"
"""
 
## External Mods
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

## User Mods
import ocean_irradiance_module.Ocean_Irradiance as OI
#import ocean_irradiance_baird.ocean_irradiance_baird as OIB
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
#import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC
from ocean_irradiance_module import Wavelength_To_RGB
import reflectance_spectra
from ocean_irradiance_module.Phytoplankton_Colormap import Get_Phy_Cmap_Dict


def Plot_Chla_vs_Rrs_Ratio(Rrs_field_species, wavelengths, species, chlas, method='shubha'):
    """
    This is to create the plot as suggested by Shubha as described in the file dcostring above.
    """
    cmap = Get_Phy_Cmap_Dict()

    fig, ax = plt.subplots()
    xlim = [min(chlas), max(chlas)]
    ## Looping the species.
    for k,phy_type in enumerate(species):
        ## Making the rrs ratio of rrs443/rrs551.
        Rrs_arr = Rrs_field_species[phy_type]
        ## 443 should be index 0, 551 index 1.
        Rrs_ratio = np.maximum(Rrs_arr[:, 0], Rrs_arr[:,1], Rrs_arr[:,2]) / Rrs_arr[:,3]
        
        ## Now for the plotting.
        ax.semilogy(Rrs_ratio, chlas, color=cmap[phy_type], label=phy_type)
#        ax.loglog(chlas, Rrs_ratio, color=colors[k], label=phy_type)
    
    ## Plotting the OCx algorithim chla.
    OCx_Rrs_ratio = np.linspace(0, 12, 100)
    ## Setting the blue to be the ratio and Rrs for green to be 1, s.t. ratio is Rrsb. 
    OC4_chla = OIR.OCx_alg(OCx_Rrs_ratio, 0, 0, np.ones(len(OCx_Rrs_ratio)), method ='OC4')
    OC3V_chla = OIR.OCx_alg(OCx_Rrs_ratio, 0, 0, np.ones(len(OCx_Rrs_ratio)), method ='OC3V')
    OC4E_chla = OIR.OCx_alg(OCx_Rrs_ratio, 0, 0, np.ones(len(OCx_Rrs_ratio)), method ='OC4E')
    OC3M_chla = OIR.OCx_alg(OCx_Rrs_ratio, 0, 0, np.ones(len(OCx_Rrs_ratio)), method ='OC3M')
    ax.semilogy(OCx_Rrs_ratio, OC4_chla, '--', color='k', label='OC4')
    ax.semilogy(OCx_Rrs_ratio, OC3V_chla, ':', color='k', label='OC3V')
    ax.semilogy(OCx_Rrs_ratio, OC4E_chla, '-.', color='k', label='OC4E')
    ax.semilogy(OCx_Rrs_ratio, OC3M_chla, linestyle=(0, (3,5,1,5,1,5)), color='k', label='OC3M')
#    ax.loglog(OCx_chla, OCx_Rrs_ratio, '--', color='k', label='OCx')
    
    ax.legend(title='Phy Species')
    ax.grid()
    ax.set_ylim(xlim)
    ax.set_xlim([0,12])
    ax.set_xlabel(r'$\frac{\mathrm{R_{rs}}(\lambda_{\mathrm{blue}})}{\mathrm{R_{rs}}(\lambda_{\mathrm{green}})}$', fontsize=14) 
    ax.set_ylabel(r'Chla [mg $\mathrm{m}^{-3}$]')
    if method == 'shoot_up':
        ax.set_title('UCSC Version of Dutkiewicz et al. (2015) \n Radiative Transfer Model with Uniform Chl-a Profiles')
    if method  == 'shubha': 
        ax.set_title('UCSC Version of Sathyendranath et al. (1997) \n Two Stream Radiative Transfer Model with Uniform Chl-a Profiles')

    fig.show()

    return 

     
if __name__ == '__main__': 
    
    ## Params 
    PI = Param_Init()
    
    ## The number of chla values 
    N_chla = 100
    ## The number of vertical levels.
    N = 200

    ## The chla values 
    chlas = np.logspace(-1, 1.5, N_chla)
#    chlas = np.linspace(0.1, 30, N_chla)
    ## The chla full array
    chla_array = np.zeros((N, len(chlas)))
    for k, chla in enumerate(chlas):
        chla_array[:,k] = np.full(N,chla)
    ## The accompanying z-grid
    z = np.linspace(-1000, 0, N)

    ## The wavelengths for the ratio 
    wavelengths = [443, 486, 510, 551]
    
    ## The phytoplankton species
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']

    ## The irradiance method to use. 
    method = 'scipy'
    
    ## Using the solve irradiance funtion from the file  reflectance_spectra.py
    Rrs_field_species = reflectance_spectra.Solve_Irradiance(PI, N,  wavelengths, species, chla_array, z, method=method)

    ## plotting the restult.
    Plot_Chla_vs_Rrs_Ratio(Rrs_field_species, wavelengths, species, chlas, method=method)




