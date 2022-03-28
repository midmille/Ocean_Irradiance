"""
Created: March 9, 2022 7:48am
Author: Miles D. Miller, University of California Santa Cruz

This file is for the purpose of running the irradaince light model using 
OOI data as well as comparing the model output to in situ data.

To run: 
    --> Be sure to be in the ooi environment. 
        -https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python
"""

## [Import Statements]
## [User Modules]
import ooi_data_functions as ODF 
from ocean_irradiance_module import Ocean_Irradiance as OI
import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module import Wavelength_To_RGB

## [External Modules]
import numpy as np
import pickle
from ooi_data_explorations.data_request import data_request

def Download_Data(site_name, assembly, instrument, method, start, stop):
    """
    This file is to download the ooi data from source using the ooi_data_explorations.data_request module. 
    
    Parameters
    ----------
    ooi_paramdict: Dict
        This is a dictionary containing the required arguments for the OOI data download. 
        The keys are as follows: [site, assembly, instrument, method, start, stop]. 
    Returns
    -------
    ooi_dat: 

    """

    print('HERE')
    print(site_name)
    ooi_dat = data_request(site_name, assembly, instrument, method, start=start, stop=stop)

    return ooi_dat


def Get_Chla_Profiles(ooi_dat):
    """
    This function retrieves list of profiles for depth, datetime, and chla. 

    Parameters
    ----------

    Returns
    -------
    """

    ## [First get the depth, dt, and chla arrays from ooi_dat.]
    depth_dat = ooi_dat.variables['depth']
    dt_dat = ooi_dat.variables['time']
    chla_dat = ooi_dat.variables['fluorometric_chlorophyll_a']

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, chla_profs = ODF.Create_Profiles(depth_dat, dt_dat, chla_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, chla_profs


def Get_Spkir_Profiles(ooi_dat):
    """
    This function retrieves list of profiles for depth, datetime, and chla. 

    Parameters
    ----------

    Returns
    -------
    """

    ## [First get the depth, dt, and chla arrays from ooi_dat.]
    depth_dat = ooi_dat.variables['depth']
    dt_dat = ooi_dat.variables['time']
    spkir_dat = ooi_dat.variables['spkir_abj_cspp_downwelling_vector']

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, spkir_profs = ODF.Create_Profiles(depth_dat, dt_dat, spkir_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, spkir_profs



def Run_Irradiance(N, wavelengths, phy_type, depth_profs, dt_profs, chla_profs):
    """
    This function runs the irradiance model over all the different chla profiles.

    """
    ## [The irradiance parameters object.]
    PI = Param_Init()

    ## [The number of profiles.]
    N_profs = len(depth_profs)

    ## [Irradiance field dictionary.]
    irr_field = {}

    ## [The irradiance array.]
    irr_arr = np.zeros((N, N_profs, 4))

    ## [Loop over the wavelengths.]
    for lam in wavelengths:
        ## [Loop over the OOI profiles.]
        for k in range(N_profs):
            ## [Make the depth negative.]
            depth = -depth_profs[k]
            ## [Make the depth profile start at zero.]
            depth = depth
            ## [Create the phytoplankton object.]
            phy = OI.Phy(depth_profs[k], chla_profs[k], esd(phy_type), 
                         abscat(lam, phy_type, C2chla='default')[0], 
                         abscat(lam, phy_type, C2chla='default')[1])
            ocean_irr_sol = OIS.ocean_irradiance_shubha(depth[0], 
                                                        PI.Ed0+PI.Es0, 
                                                        abscat(lam, 'water'), 
                                                        PI.coefficients, 
                                                        phy=phy, 
                                                        CDOM=None, 
                                                        N=N, 
                                                        pt1_perc_zbot=False, 
                                                        pt1_perc_phy=False)
            ## [Storing output into array.]
            irr_arr[:,k,0] = ocean_irr_sol[0]
            irr_arr[:,k,2] = ocean_irr_sol[1]
            irr_arr[:,k,3] = ocean_irr_sol[2] 
        
        ## [Store array to dictionary.]
        irr_field[lam] = irr_arr

    return irr_field


def Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_wavelengths, irr_field, spkir_depth_profs, spkir_dt_profs, spkir_profs): 
    """
    This function plots downwelling irradiance stream solved for using the chla profiles
    from ooi. Then it plots the spkir data that corresponds to the solved irradiance profiles.

    Parameters
    ----------
    prof_index: Integer
        The index of the profile to be plotted. 
    wavelengths: List, Integer
        The list of wavelengths, the wavelengths of the ooi SPKIR might differ from the
        wavelengths used for the coefficients in the irradiance model. Might has to pass
        the spkir wavelengths eventually. 
    irr_field: Dictionary, 4-D Array
        The irradiance field dictionary with wavelengths as keys and the 4-D irradiance array
        as the values.
    spkir_depth_profs: List, 1-D Array
        The list of spkir depth profiles. The depth profiles are positive with the bottom most
        depth at index 0 and the surface at index -1. 
    spkir_dt_profs: List, 1-D Array
        The list of the times associated with each depth and spkir data point. 
    spkir_profs: List, 2-D Array
        The list of spkir profiles. The first index of each lsit value is the profile and the
        second index corresponds to the wavelenth.
    """

    ## [The number of wavelengths.]    
    N_lam = len(wavelengths)

    ## [The spkir data for the given index.]
    depth = -spkir_depth_profs[prof_index]
    dt = spkir_dt_profs[prof_index]
    spkir = spkir_profs[prof_index]

    fig, ax = plt.subplots()
    ## [Get the colors such that they match the wavelengths.] 
    colors = [Wavelength_To_RGB.wavelength_to_rgb(wavelength) for wavelength in wavelengths]
    
    ## [Loop over the irradiance wavelengths.]
    for k, lam in enumerate(wavelengths): 
        irr_arr = irr_field[lam]
        ## [ Plotting the downward direct profile. 0 is downward irradiance, 3 is depth.]
        ax.plot(irr_arr[:, prof_index, 0], irr_arr[:, prof_index, 0], ':', label=f'Irr {lam}', color=colors[k])

    for k, lam in enumerate(spkir_wavelengths):
        ## [Plotting the spkir profile.]
        
    

    

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Run irradiance model using OOI profile input.')
#    parser.add_argument('--ooi_paramdict', help='A dictionary containing the necessary arguments to'\
#                                               ' download the ooi data. The keys to the dictionary'\
#                                               ' should be the following: [site, assembly, instrument,'\
#                                               ' method, start, stop]')
    parser.add_argument('--site_name', help='site name')
    parser.add_argument('--assembly', help= 'assembly')
    parser.add_argument('--method', help = 'method')
    parser.add_argument('--start', help='start')
    parser.add_argument('--stop', help='stop')
    

    
    parser.add_argument('--ooi_savefile_head', help='The name of the pickle file in which the ooi data is saved.'\
                                               'If this argument is given then the ooi data will be loaded'\
                                               ' from file and not downloaded from source.', 
                                               type=str, required=True)
    args = parser.parse_args()

    ## [The savefile names for each ooi data set.]
    spkir_savefile = args.ooi_savefile_head + '_spkir.p'
    flort_savefile = args.ooi_savefile_head + '_flort.p'

    ## [Download OOI data.]
    if args.site_name: 
        ## [Download the data.]
        ## [The Chla data.]
        flort_dat = Download_Data(args.site_name, args.assembly, 'FLORT', args.method, args.start, args.stop)
        ## [The spkir data.]
        spkir_dat = Download_Data(args.site_name, args.assembly, 'SPKIR', args.method, args.start, args.stop)
        ## [Save the ooi data into the pickle file.]
        pickle.dump(flort_dat, open(flort_savefile, 'wb'))
        pickle.dump(spkir_dat, open(spkir_savefile, 'wb'))

    ## [Load the data from pickle.]
    else:
        flort_dat = pickle.load(open(flort_savefile, 'rb'))
        spkir_dat = pickle.load(open(spkir_savefile, 'rb'))

    ## [Get the chla profile lists.]
    depth_profs, dt_profs, chla_profs = Get_Chla_Profiles(flort_dat)
    ## [Get the spkir profile lists.]
    spkir_depth_profs, spkir_dt_profs, spkir_profs = Get_Spkir_Profiles(spkir_dat)

    ## [Get the spkir wavelengths.]
    spkir_wavelengths = ODF.Get_SPKIR_Wavelengths(spkir_dat.variables['spkir_abj_cspp_downwelling_vector'])

    ## [Run the irradiance model using the profiles, over all profiles.]
    N=100
    wavelengths = [443]
    phy_type = 'Diat'
    irr_field = Run_Irradiance(N, wavelengths, phy_type, depth_profs, dt_profs, chla_profs)

