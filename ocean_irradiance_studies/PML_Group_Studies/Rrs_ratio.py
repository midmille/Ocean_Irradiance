"""
Created: January 31 2021
Author: Miles D Miller 


This file follows Shuhba's inquiry from an email on January 21, 2022, 5:56am

"It would also be good to plot on a log-log scale chl on the x-axis and 
Rrs ratio (say Rrs443/Rrs555) or the closest wavelength pair that you have,
again for individual phytoplankton types"
"""
 
## User Mods
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_baird.ocean_irradiance_baird as OIB
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC
from ocean_irradiance_module import Wavelength_To_RGB
import reflectance_spectra


def


if __name__ == '__main__': 
    
    ## The number of chla values 
    N_chla = 300
    ## The number of vertical levels.
    N = 200

    ## The chla values 
    chlas = np.linspace(.01, 10, N_chla)
    ## The chla full array
    chla_array = np.zeros((N, len(chlas)))
    ## Making a full profile for each chla value. 
    chla_array = np.zeros((N, len(chlas)))
    for k, chla in enumerate(chlas):
        chla_array[:,k] = np.full(N,chla)
    ## The accompanying z-grid
    z = np.linspace(-1000, 0, N)

    ## The wavelengths for the ratio 
    wavelengths = [443, 551]
    
    ## The phytoplankton species
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']
    
    ## The pickle file
    #save_file = 

