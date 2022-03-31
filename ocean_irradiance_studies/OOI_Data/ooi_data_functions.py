"""
Created: 10:59am February 25th 2022
Author: Miles D. Miller

This file contains some useful fucntions for parsing over the ooi data. 
"""

## [External Mods]
from ooi_data_explorations.data_request import data_request
import matplotlib.pyplot as plt
import re
import numpy as np
## [User Mods]
from ocean_irradiance_module import Wavelength_To_RGB


def Create_Profiles(depth, date_time, val, dz_max): 
    """
    This seperates the 1-D Array for the depth and corresponding value 
    into a list of Smaller 1-D Arrays which are just a single profile

    Parameters
    ----------
    depth: 1-D Array [N]
        The depths corrspnding to the given values. 
    val: 1-D Array [N]
        The values of the data at each depth to be split into a list
        of 1-D profiles. 
    dz_max: Float
        The maximum value of separation to decide where one profile begins
        and another ends. 

    Returns
    -------
    depth_profs: List
        The list of depth arrays for each profile. 
    val_profs: List
        The list of values for each profile. 
    """

    ## [The start of the function.]
    ## [Init the list.]
    depth_profs = []
    date_time_profs  = []
    val_profs = []
    ## [Init the index that marks where the profile begins.]
    k_prof = 0

    ## [Loop over the entire 1-D depth array.]
    for k in range(len(depth)): 
    
        ## [Check if not on the last k.]
        if (k != len(depth) - 1): 

            ## [This statement determines if at the end of profile or not.]
            if abs(depth[k+1] - depth[k]) > dz_max: 
                ## [The array of the profile.]
                depth_prof = depth[k_prof:k]
                date_time_prof = date_time[k_prof:k]
                val_prof = val[k_prof:k]
                ## [Appending the arrays to their respective lists.]
                depth_profs.append(depth_prof)
                date_time_profs.append(date_time_prof)
                val_profs.append(val_prof)
                ## [Rewriting the index where the profiles ends/ next begins.]
                k_prof = k+1

        ## [If loop is at last index.]
        if (k == len(depth)-1): 
            ## [The array of the profile.]
             depth_prof = depth[k_prof:k]
             date_time_prof = date_time[k_prof:k]
             val_prof = val[k_prof:k]
             ## [Appending the arrays to their respective lists.]
             depth_profs.append(depth_prof)
             date_time_profs.append(date_time_prof)
             val_profs.append(val_prof)
            
                
    return depth_profs, date_time_profs, val_profs   
            

def Get_SPKIR_Wavelengths(spkir_dat):
    """
    This gets the wavelengths used in the ooi spkir data set from the attributes comment

    Parameters
    ----------
    spkir_dat: ooi data_request object
        The spkir data set. 

    Returns
    -------
    spkir_wavelengths: List, Integer
        The list of the wavelengths from the given SPKIR data set. 
    """

    ## [Get the wavelengths from the comment for the spkir data.]
    spkir_comment =  spkir_dat.attrs['comment']
    ## [Find the locations of the 'nm' units in the comment.]
    wavelength_i = [_.start() for _ in re.finditer('nm', spkir_comment)] 
    ## [The wavelengths list init.]
    spkir_wavelengths = []
    ## [Loop over the 'nm' instances adn get the wavelengths']
    for k in wavelength_i: 
        spkir_wavelengths.append(int(spkir_comment[k-3:k]))
        
    return spkir_wavelengths


def Plot_SPKIR_Profile(prof_index, spkir_data, site, assembly, instrument, method): 
    """

    """

    ## [Read in the diffferent variables for spkir data set.]
    depth_dat = spkir_data.variables['depth']
    date_time_dat = spkir_data.variables['time']
    spkir_dat = spkir_data.variables['spkir_abj_cspp_downwelling_vector']

    ## [Get the corresponding data in profile list form.]
    dz_max = 1
    depth_profs, date_time_profs, spkir_profs = Create_Profiles(depth_dat.data, date_time_dat.data, spkir_dat.data, dz_max)
    ## [Make the given profs theri own variables.]
    depth = depth_profs[prof_index]
    for d in depth:
        print(d)
    date_time = date_time_profs[prof_index]
    spkir = spkir_profs[prof_index]

    wavelengths = Get_SPKIR_Wavelengths(spkir_dat)
    
    ## [Plot the profile for th egiven prof_index.]
    fig, ax = plt.subplots()
    ## [Plot the profile for each wavelength.]
    colors = [Wavelength_To_RGB.wavelength_to_rgb(wavelength) for wavelength in wavelengths]
    for k in range(len(wavelengths)):
        ax.plot(spkir[:,k], -depth, label = wavelengths[k], color = colors[k])
    ## [Labels.]
    ax.set_ylabel(f"Z [{depth_dat.attrs['units']}]")
    ax.set_xlabel(f"Downwelling Spectral Irradiance {spkir_dat.attrs['units']}")
    ax.set_title(f"OOI SPKIR Profile from \n {date_time[0]} to {date_time[-1]}")
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
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: {instrument}')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.legend(title='Wavelengths [nm]')
    ax.grid()

    fig.show()

    return 


def Plot_OOI_Profile(prof_index, data_set, var_name, site, assembly, instrument, method): 
    """

    """

    ## [Read in the diffferent variables for spkir data set.]
    depth_dat = data_set.variables['depth']
    dt_dat = data_set.variables['time']
    var_dat = data_set.variables[var_name]

    ## [Get the corresponding data in profile list form.]
    dz_max = 1
    depth_profs, dt_profs, var_profs = Create_Profiles(depth_dat.data, dt_dat.data, var_dat.data, dz_max)
    ## [Make the given profs theri own variables.]
    depth = -depth_profs[prof_index]
    dt = dt_profs[prof_index]
    var = var_profs[prof_index]

    ## [Plot the profile for th egiven prof_index.]
    fig, ax = plt.subplots()
    ax.plot(var, depth)
    ## [Labels.]
    ax.set_ylabel(f"Z [{depth_dat.attrs['units']}]")
    ax.set_xlabel(f"{var_name} {var_dat.attrs['units']}")
    ax.set_title(f"OOI {var_name} from \n {dt[0]} to {dt[-1]}")
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
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: {instrument}')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.grid()

    fig.show()

def Plot_Chla_Hoffmuller(flort_data, site, assembly, instrument, method): 
    """

    """

    ## [Read in the diffferent variables for spkir data set.]
    depth_dat = flort_data.variables['depth']
    date_time_dat = flort_data.variables['time']
    chla_dat = flort_data.variables['fluorometric_chlorophyll_a']

    ## [Get the corresponding data in profile list form.]
    dz_max = 1
    depth_profs, date_time_profs, chla_profs = Create_Profiles(depth_dat.data, date_time_dat.data, chla_dat.data, dz_max)

    ## [Create the hoffmuller diagram.]
    ## [First find the largest size profile in the list of profiles.]
    ## [This largest size will set the vertical size of the mesh.]
    sizes = np.zeros(len(chla_profs))
    for k,chla_prof in enumerate(chla_profs): 
        sizes[k] = len(chla_prof)
    size = max(sizes)
        
    ## [Constructing the mesh.]
    chla_mesh = np.NaN * np.zeros((int(size), len(chla_profs)))
    ## [The start times of each prof for x-coords.]
    t_coord = np.zeros(len(date_time_profs))

    ## [Filling the mesh with profiles, notes they are note interpolated
    ##  so this method is incorrect, but it will work for now.]
    for k, chla_prof in enumerate(chla_profs):
        prof_len = len(chla_prof)
        chla_mesh[:prof_len, k] = chla_prof
        t_coord[k] = date_time_profs[k][0] 
        if len(chla_prof) == size: 
            z_coord =  depth_profs[k]

    ## [Plotting the Hovmuller disgram.]
    fig, ax=plt.subplots()
    tt, zz = np.meshgrid(t_coord, z_coord)
    im = ax.pcolormesh(tt, zz, chla_mesh, cmap='nipy_spectral')
    fig.colorbar(im, ax=ax)

    fig.show()

    return 






