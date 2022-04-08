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
import matplotlib.pyplot as plt

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


def Get_Flort_Profiles(ooi_dat):
    """
    This function retrieves list of profiles for depth, datetime, chla, and opticalbackcatter.

    Parameters
    ----------

    Returns
    -------
    """

    ## [First get the depth, dt, and chla arrays from ooi_dat.]
    depth_dat = ooi_dat.variables['depth']
    dt_dat = ooi_dat.variables['time']
    chla_dat = ooi_dat.variables['fluorometric_chlorophyll_a']
    opt_backscatter_dat = ooi_dat.variables['optical_backscatter']

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, chla_profs = ODF.Create_Profiles(depth_dat, dt_dat, chla_dat, dz_max)    
    depth_profs, dt_profs, opt_bs_profs = ODF.Create_Profiles(depth_dat, dt_dat, opt_backscatter_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, chla_profs, opt_bs_profs


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


def Get_Optaa_Profiles(optaa_dat): 
    """
    This function retireves the profiles OPTAA absorption data.
    """
    ## [First get the depth, dt, and chla arrays from ooi_dat.]
    depth_dat = optaa_dat.variables['depth']
    dt_dat = optaa_dat.variables['time']
    abs_dat = optaa_dat.variables['optical_absorption']
    wavelength_dat = optaa_dat.variables['wavelength_a']

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, abs_profs = ODF.Create_Profiles(depth_dat, dt_dat, abs_dat, dz_max)    
    depth_profs, dt_profs, wavelength_profs = ODF.Create_Profiles(depth_dat, dt_dat, wavelength_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, abs_profs, wavelength_profs
 



def Run_Irradiance(N, wavelengths, phy_type, depth_profs, dt_profs, chla_profs, opt_bs_profs, optaa_dat):
    """
    This function runs the irradiance model over all the different chla profiles.

    """
    ## [The irradiance parameters object.]
    PI = Param_Init()

    ## [The number of profiles.]
    #N_profs = len(depth_profs)
    N_profs = 1

    ## [Irradiance field dictionary.]
    irr_field = {}
    irr_field_ab = {}

    ## [Loop over the wavelengths.]
    for lam in wavelengths:
        ## [The irradiance array.]
        irr_arr = np.zeros((N, N_profs, 4))
        irr_arr_ab = np.zeros((N, N_profs, 4))

        ## [Loop over the OOI profiles.]
        for k in range(N_profs):
            ## [Make the depth negative.]
            depth = -depth_profs[k]
            ## [Make the depth profile start at zero.]
            print('Chla_prof[k]', chla_profs[k])
            ## [Create the phytoplankton object.]
            phy = OI.Phy(depth, chla_profs[k], esd(phy_type), 
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
            a_irr = ocean_irr_sol[3]
            b_b_irr = ocean_irr_sol[4]

            print('a_irr:', a_irr)
            print('b_b_irr:', b_b_irr)


            ## [The optaa data.]
            optaa_depth_dat = optaa_dat.variables['depth'].data
            optaa_dt_dat = optaa_dat.variables['time'].data
            optaa_abs_dat = optaa_dat.variables['optical_absorption'].data
            optaa_wavelength_dat = optaa_dat.variables['wavelength_a'].data

            dt_lbnd = dt_profs[k].data[0] 
            dt_ubnd = dt_profs[k].data[-1]
            prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, optaa_dt_dat)

            optaa_dt_prof = optaa_dt_dat[prof_mask]
            optaa_depth_prof = optaa_depth_dat[prof_mask]
            optaa_abs_prof = optaa_abs_dat[prof_mask]
            optaa_wavelengths = optaa_wavelength_dat[prof_mask]

            print(optaa_wavelengths.shape)

            optaa_depth = -np.squeeze(optaa_depth_prof)
            ## [Must get the wavelength index.]
            lam_i = ODF.Get_Wavelength_Index(optaa_wavelengths[0,:], lam)
            print('lam_i',lam_i)
            ## [The squeese makes the array 1-D.]
            a = np.squeeze(optaa_abs_prof[:, lam_i])
            b = opt_bs_profs[k]

            
            ## [Plotting the absorption comp]
            fig, ax = plt.subplots()
            ax.plot(a, optaa_depth, label='Abs OOI OPTAA')
            ax.plot(a_irr, irr_arr[:,k,3], label='Abs Irr')
            ax.set_title('Absorption')
            ax.set_ylabel('Z[m]')
            ax.set_xlabel('Abs [m^-1]')
            ax.legend()
            ax.grid()
            fig.show()

            ## [Plotting the scattering comp.]
            fig,ax = plt.subplots()
            ax.plot(b, depth, label='Scat OOI FLORT')
            ax.plot(b_b_irr, irr_arr[:,k,3], label='Scat Irr')
            ax.set_title('Back Scattering')
            ax.set_ylabel('Z[m]')
            ax.set_xlabel('Scat [m^-1]')
            ax.legend()
            ax.grid()
            fig.show()
            


            print(a.shape)

            ## [Now running the irr model using abs and scat from ooi.]
            ocean_irr_sol_ab = OIS.ocean_irradiance_shubha_ab(depth[0], 
                                                              PI.Ed0+PI.Es0, 
                                                              PI.coefficients, 
                                                              optaa_depth.data,  ## [The z_a is optta depth.]
                                                              a.data, 
                                                              depth.data, ## [The z_b is the flort depth.]
                                                              b.data,
                                                              CDOM=None, 
                                                              N=N, 
                                                              pt1_perc_zbot=False)


            ## [Storing output into array.]
            irr_arr_ab[:,k,0] = ocean_irr_sol_ab[0]
            irr_arr_ab[:,k,2] = ocean_irr_sol_ab[1]
            irr_arr_ab[:,k,3] = ocean_irr_sol_ab[2] 
        
        ## [Store array to dictionary.]
        irr_field[lam] = irr_arr
        irr_field_ab[lam] = irr_arr_ab

    return irr_field, irr_field_ab


def Irr_Abs_Scat(): 
    """
    This function returns the absorpostion and scttering as would be caclulated by the irradiance function.
    """

    return 

def Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_wavelengths, spkir_wavelength_index,irr_field, irr_field_ab, spkir_dat, site, assembly, method): 

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

    ## [spkir depth stuff.]
    spkir_depth_dat = spkir_dat.variables['depth'].data
    spkir_dt_dat = spkir_dat.variables['time'].data
    spkir_dat = spkir_dat.variables['spkir_abj_cspp_downwelling_vector'].data

    dt_lbnd = dt_profs[prof_index].data[0] 
    dt_ubnd = dt_profs[prof_index].data[-1]
    prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, spkir_dt_dat)

    spkir_dt_prof = spkir_dt_dat[prof_mask]
    spkir_depth_prof = spkir_depth_dat[prof_mask]
    spkir = spkir_dat[prof_mask]

    depth = -spkir_depth_prof 

    fig, ax = plt.subplots()
    ## [Get the colors such that they match the wavelengths.] 
    colors = [Wavelength_To_RGB.wavelength_to_rgb(wavelength) for wavelength in wavelengths]
    
    ## [Loop over the irradiance wavelengths.]
    for k, lam in enumerate(wavelengths): 
        irr_arr = irr_field[lam]
        irr_arr_ab = irr_field_ab[lam]
        ## [ Plotting the downward direct profile. 0 is downward irradiance, 3 is depth.]
        ## [Also must multiply by surface value of spkir prof, since irr_arr is normalized.]
        ax.plot(spkir[:,spkir_wavelength_index[k]].data[-1] * irr_arr[:, prof_index, 0], irr_arr[:, prof_index, 3], ':', label=f'Irr {lam}', color=colors[k])
        ax.plot(spkir[:,spkir_wavelength_index[k]].data[-1] * irr_arr_ab[:, prof_index, 0], irr_arr_ab[:, prof_index, 3], '-', label=f'Irr ooi ab {lam}', color=colors[k])

    for k,i in enumerate(spkir_wavelength_index):
        ## [Plotting the spkir profile.]
        ## [Make the 1m avg grid.]
        depth_avg = np.arange(depth[0], depth[-1], 1)
        spkir_avg = ODF.Grid_Average_Profile(depth, spkir[:,i], depth_avg)
        #ax.plot(spkir[:, i], depth, '--', label=f'OOI SPKIR {lam}', color=colors[k])
        ax.plot(spkir_avg, depth_avg, '--', label=f'OOI SPKIR {lam}', color=colors[k])

    ## [Labels.]
    #ax.set_ylabel(f"Z [{depth_dat.attrs['units']}]")
    ax.set_ylabel(f"Z [m]")
    #ax.set_xlabel(f"Downwelling Spectral Irradiance {spkir_dat.attrs['units']}")
    ax.set_xlabel(f"Downwelling Spectral Irradiance")
    #ax.set_title(f"OOI SPKIR Profile and Irradiance Model \n OOI SPKIR Profile Date: {date_time[0]} to {date_time[-1]}")
    ax.set_title(f"OOI SPKIR Profile and Irradiance Model")
    ## [Putting some identifying text on the figure.]
    ## [10% up the vertical location]
    txt_y = ax.get_ylim()[1] + 0.5 * ax.get_ylim()[0] 
    ## [10% of the horizontal location.]
    txt_x = ax.get_xlim()[0] + 0.2 * ax.get_xlim()[1]
    ## [The change in txt location in vertical.]
    txt_dz = 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ## [Adding the txt.]
    ax.text(txt_x, txt_y, f'SITE: {site}')   
    ax.text(txt_x, txt_y+txt_dz, f'ASSEMBLY: {assembly}')   
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: SPKIR')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.legend(title='Wavelengths [nm]')
    ax.grid()

    fig.show()


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Run irradiance model using OOI profile input.')
#    parser.add_argument('--ooi_paramdict', help='A dictionary containing the necessary arguments to'\
#                                               ' download the ooi data. The keys to the dictionary'\
#                                               ' should be the following: [site, assembly, instrument,'\
#                                               ' method, start, stop]')
    parser.add_argument('--download', action='store_true', help='Download the data from OOI [This takes time], must specify file head')
    parser.add_argument('--load', action='store_true', help='Load the data from pickle, must specify file head.')
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
    optaa_savefile = args.ooi_savefile_head + '_optaa.p'

    ## [Download OOI data.]
    if args.download: 
        ## [Download the data.]
        ## [The Chla data.]
        flort_dat = Download_Data(args.site_name, args.assembly, 'FLORT', args.method, args.start, args.stop)
        ## [The spkir data.]
        spkir_dat = Download_Data(args.site_name, args.assembly, 'SPKIR', args.method, args.start, args.stop)
        ## [The OPTAA data.]
        optaa_dat = Download_Data(args.site_name, args.assembly, 'OPTAA', args.method, args.start, args.stop)
        ## [Save the ooi data into the pickle file.]
        pickle.dump(flort_dat, open(flort_savefile, 'wb'))
        pickle.dump(spkir_dat, open(spkir_savefile, 'wb'))
        pickle.dump(optaa_dat, open(optaa_savefile, 'wb'))

    ## [Load the data from pickle.]
    elif args.load:
        flort_dat = pickle.load(open(flort_savefile, 'rb'))
        spkir_dat = pickle.load(open(spkir_savefile, 'rb'))
        optaa_dat = pickle.load(open(optaa_savefile, 'rb'))

    ## [Get the chla profile lists.]
    depth_profs, dt_profs, chla_profs, opt_bs_profs = Get_Flort_Profiles(flort_dat)

    ## [Get the spkir wavelengths.]
    spkir_wavelengths = np.array(ODF.Get_SPKIR_Wavelengths(spkir_dat.variables['spkir_abj_cspp_downwelling_vector']))

    ## [Run the irradiance model using the profiles, over all profiles.]
    N=100
    wavelengths = [443]
    phy_type = 'Generic'
    irr_field, irr_field_ab = Run_Irradiance(N, wavelengths, phy_type, depth_profs, dt_profs, chla_profs, opt_bs_profs, optaa_dat)

    ## [Some params for plotting.]
    prof_index = 0 
    ## [The index of the spkir wavelengths that most closely matches the given irradiance wavcelengths.]
    ## [This wavelength index should be automated soon.]
    spkir_wavelength_index = [ODF.Get_Wavelength_Index(spkir_wavelengths, wavelengths[0])]
    print('lam_i',spkir_wavelength_index)
    ## [Plot the OOI SPKIR against the irradiance model.]
    Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_wavelengths, spkir_wavelength_index, irr_field, irr_field_ab, spkir_dat, args.site_name, args.assembly, args.method)


