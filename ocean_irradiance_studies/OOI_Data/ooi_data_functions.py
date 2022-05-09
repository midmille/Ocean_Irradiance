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
import pickle
## [User Mods]
from ocean_irradiance_module import Wavelength_To_RGB


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

    ooi_dat = data_request(site_name, assembly, instrument, method, start=start, stop=stop)

    return ooi_dat


def Download_OOI_Data(ooi_savefile_head, download, site_name, assembly, method, start, stop):
    """
    This function is for the download xor load of the desired OOI data. If the files already exists then pass
    the download argument as false. Otherwise if you want to download the data from OOI servers, 
    which takes some times, then pass the download argument as true.

    Parameters
    ----------
    ooi_savefile_head: String
        The header of the save file. The function will add '_{instrumentname}.p' for each different instrument, 
        which is [spkir, flort, optaa]. 
    download: Boolean 
        Set True to download OOI data from servers. Set False to load from predownloaded pickle files.
    site_name: String
        The site name variables for the data.
    assembly: String 
        The assembly. 
    method: String 
        The method.
    start: String
        The start date of the profiler data. 
    stop: String
        The stop date of the profiler data. 

    Returns
    -------
    flort_dat: ooi data_request object
        The flort data.
    spkir_dat: ooi data_request object
        The spkir data.
    optaa_dat: ooi data_request object
        The optaa data.

    """
    ## [The savefile names for each ooi data set.]
    spkir_savefile = ooi_savefile_head + '_spkir.p'
    flort_savefile = ooi_savefile_head + '_flort.p'
    optaa_savefile = ooi_savefile_head + '_optaa.p'

    ## [Download OOI data.]
    if download: 
        ## [Download the data.]
        ## [The Chla data.]
        flort_dat = Download_Data(site_name, assembly, 'FLORT', method, start, stop)
        ## [The spkir data.]
        spkir_dat = Download_Data(site_name, assembly, 'SPKIR', method, start, stop)
        ## [The OPTAA data.]
        optaa_dat = Download_Data(site_name, assembly, 'OPTAA', method, start, stop)
        ## [Save the ooi data into the pickle file.]
        pickle.dump(flort_dat, open(flort_savefile, 'wb'))
        pickle.dump(spkir_dat, open(spkir_savefile, 'wb'))
        pickle.dump(optaa_dat, open(optaa_savefile, 'wb'))

    ## [Load the data from pickle.]
    else:
        flort_dat = pickle.load(open(flort_savefile, 'rb'))
        spkir_dat = pickle.load(open(spkir_savefile, 'rb'))
        optaa_dat = pickle.load(open(optaa_savefile, 'rb'))

    return flort_dat, spkir_dat, optaa_dat


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


def Get_Flort_Dat(flort_dat): 
    """
    This retrieves all the useful flort data. 
    """

    ## [Make the depth negative.]
    depth_dat = - flort_dat.variables['depth']
    dt_dat = flort_dat.variables['time']
    chla_dat = flort_dat.variables['fluorometric_chlorophyll_a']
    cdom_dat = flort_dat.variables['fluorometric_cdom']
    opt_bs_dat = flort_dat.variables['optical_backscatter']

    return depth_dat, dt_dat, chla_dat, cdom_dat, opt_bs_dat


def Get_Flort_Profiles(flort_dat):
    """
    This function retrieves list of profiles for depth, datetime, chla, and opticalbackcatter.

    Parameters
    ----------

    Returns
    -------
    """

    ## [Retrieve the data.]
    depth_dat, dt_dat, chla_dat, cdom_dat, opt_bs_dat = Get_Flort_Dat(flort_dat)

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, chla_profs = Create_Profiles(depth_dat, dt_dat, chla_dat, dz_max)    
    depth_profs, dt_profs, cdom_profs = Create_Profiles(depth_dat, dt_dat, cdom_dat, dz_max)    
    depth_profs, dt_profs, opt_bs_profs = Create_Profiles(depth_dat, dt_dat, opt_bs_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs


def Get_Spkir_Data(spkir_dat):
    """
    This retrives all the useful SPKIR data.
    """
    
    ## [This gets all the useful data.]
    ## [Depth is made negative.]
    depth_dat = - spkir_dat.variables['depth']
    dt_dat = spkir_dat.variables['time']
    spkir_dat = spkir_dat.variables['spkir_abj_cspp_downwelling_vector']

    return depth_dat, dt_dat, spkir_dat


def Get_Spkir_Profiles(spkir_dat):
    """
    This function retrieves list of profiles for depth, datetime, and chla. 

    Parameters
    ----------

    Returns
    -------
    """

    depth_dat, dt_dat, spkir_dat = Get_Spkir_Data(spkir_dat)

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, spkir_profs = Create_Profiles(depth_dat, dt_dat, spkir_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, spkir_profs


def Get_Optaa_Dat(optaa_dat, abs_cor=False): 
    """
    This function gets the important optaa data.

    Parameters
    ---------- 
    optaa_dat: ooi_data_request data object. 
        Instrument: OPTAA, 
    abs_cor: Boolean 
        The absorption correction flag. If true the correction algorithm is implemented
        and the edited array is returned.
    """
        
    ## [The depth should be made negative.]
    depth_dat = - optaa_dat.variables['depth']
    dt_dat = optaa_dat.variables['time']
    abs_dat = optaa_dat.variables['optical_absorption']
    beam_att_dat = optaa_dat.variables['beam_attenuation']
    wavelength_dat = optaa_dat.variables['wavelength_a']
    wavelength_beam_dat = optaa_dat.variables['wavelength_c']

    ## [The empty scatter array. The underscore '_dat' is removed because the array is ]
    scat = np.zeros(abs_dat.data.shape)
    ## [Interpolate the beam attenuation onto the absorption spectral grid]
    ## [Loop over the z-coordinate grid and perform the interpolation at every z.]
    for k in range(len(depth_dat)): 
        ## [The 1-D beam attenuation array interpolated to the absroption wavelengths].
        beam_att_interp = np.interp(wavelength_dat.data[k,:], 
                                    wavelength_beam_dat.data[k,:], 
                                    beam_att_dat.data[k,:])
        ## [The total scattering is given by the total attenuation minus the absorption.]
        scat[k,:] = beam_att_interp - abs_dat.data[k,:]


    if abs_cor: 
        wavelengths = wavelength_dat.data[0,:] 
        OPTAA_Abs_Correction_G(abs_dat, wavelengths) 

    return depth_dat, dt_dat, abs_dat, scat, wavelength_dat


def Get_Optaa_Profiles(optaa_dat, abs_cor=False): 
    """
    This function retireves the profiles OPTAA absorption data.
    """

    ## [Get the optaa data.]
    depth_dat, dt_dat, abs_dat, wavelength_dat = Get_Optaa_Dat(optaa_dat, abs_cor=abs_cor)

    ## [Next retrieve profile lists using ooi_data_functions.Create_Profiles()]
    dz_max = 1
    depth_profs, dt_profs, abs_profs = ODF.Create_Profiles(depth_dat, dt_dat, abs_dat, dz_max)    
    depth_profs, dt_profs, wavelength_profs = ODF.Create_Profiles(depth_dat, dt_dat, wavelength_dat, dz_max)    

    ## [Return the profiles.]
    return depth_profs, dt_profs, abs_profs, wavelength_profs


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


def Get_Wavelength_Index(var_wavelengths, wavelength):
    """
    This function returns the wavelength index where the desired input wavelength is closest to the wavelength
    in the var_wavelengths array/list. 
    """

    diff = abs(var_wavelengths - wavelength)
    #print('diff', diff)
    wavelength_i = np.nanargmin(diff)

    return wavelength_i


def OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, var_dt): 
    """
    This function indexes the given data using the given time bounds.

    Parameters
    ----------
    
    """
    
    ## [The lower bound logic array.]
    lb = dt_lbnd < var_dt
    ## [Upper bound.]
    ub = dt_ubnd > var_dt
    ## [Total bound.]
    b = np.logical_and(lb,ub)

    return b 


def OPTAA_Abs_Correction_G(optical_absorption, wavelengths):
    """
    This function implements the correction to all absorption values in a given OPTAA
    data set. It implements the method suggested by Chris Wingard in an email on
    April 8th, 2022.

    This function implements the 715nm correction to all values of absorption in the
    data set and does not use profiles. It takes the smallest absorbtion at 715 and 
    corrects off of that value.

    Parameters
    -----------
    optical_absorption: 2-D Array, [N,M]
        optical_absorption values to be corrected. 
    wavelengths: 1-D Array, [M]
        

    Returns 
    --------
    optical_absorption: 2-D Array, [N,M]
        The optical_absorption corrected via the 715nm wavelength.
        This is the optaa data set. This is the data set with absorption corrected for
        negative values using 715nm.

    """

    ## [The 715 nm index.]
    lam_i = Get_Wavelength_Index(wavelengths, 715)

    ## [Finding all the absorption values at 715nm at all depths at all time.]
    a_715 = optical_absorption[:, lam_i]   
    
    ## [Finding the minimum at 715nm and setting that as the correction value.]
    a_cor = np.min(a_715) 

    ## [Correcting the absorption in the data set to reflect this correction.]
    optical_absorption = optical_absorption - a_cor

    return optical_absorption


def OPTAA_Abs_Correction_L(abs_profs, wavelengths):
    """
    This function implements the correction to all absorption values in a given OPTAA
    data set. It implements the method suggested by Chris Wingard in an email on
    April 8th, 2022.

    This function implements the 715nm correction to all values of absorption in the
    data set and does not use profiles. It takes the smallest absorbtion at 715 and 
    corrects off of that value.

    Parameters
    -----------
    optical_absorption: 2-D Array, [N,M]
        optical_absorption values to be corrected. 
    wavelengths: 1-D Array, [M]
        

    Returns 
    --------
    optical_absorption: 2-D Array, [N,M]
        The optical_absorption corrected via the 715nm wavelength.
        This is the optaa data set. This is the data set with absorption corrected for
        negative values using 715nm.

    """

    ## [The 715 nm index.]
    lam_i = Get_Wavelength_Index(wavelengths, 715)

    ## [Finding all the absorption values at 715nm at all depths at all time.]
    a_715 = optical_absorption[:, lam_i]   
    
    ## [Finding the minimum at 715nm and setting that as the correction value.]
    a_cor = np.min(a_715) 

    ## [Correcting the absorption in the data set to reflect this correction.]
    optical_absorption = optical_absorption - a_cor

    return optical_absorption


def Grid_Average_Profile(x_dat, y_dat, x_avg): 
    """
    This function calculates the average y-values between the cells in the x_avg grid. 
    This will help to smooth out noisy profiles. Its purpose is for plotting not for a sounds
    interpolation method. 

    x_dat and x_avg should both be increasing in value. 

    Parameters
    ----------
    x_dat: 1-D Array
        The c-coordinates of original grid the y-coordinates are on. 
    y_dat: 1-D Array
        The input data values to be smoothed out. 
    x_avg: 1-D Array
        The new grid the average values should be put on. 

    Returns
    -------
    y_avg: 1-D Array
        The average of the y-values located between cells of the x_avg grid.  
    """

    ## [Empty y_avg array.]
    ## [If the array is one or two D. ]
    if len(y_dat.shape) == 2: 
        y_avg = np.zeros((len(x_avg), y_dat.shape[1]))
        dim = 2
    elif len(y_dat.shape) == 1: 
        y_avg = np.zeros(len(x_avg))
        dim = 1
    
    if dim == 1:
        ## [Loop over the x_avg grid.]
        for k in range(len(x_avg)-1): 
            ## [y-vals when x-vals larger than the lower bound cell.]
            y_lbnd = y_dat[x_dat>x_avg[k]]
            x_lbnd = x_dat[x_dat>x_avg[k]]
            ## [y-vals smaller than the upper bound cell.]
            y_bnd = y_lbnd[x_lbnd<x_avg[k+1]]
            ## [Set the y-val avg]
            y_avg[k] = np.average(y_bnd)
        ## [Duplicate the last entry of y_avg.]
        y_avg[-1] = y_avg[-2]

    ## [if the array is 2D the routine assumes that the first dimension is the 
    ##  coordinate axis.]
    if dim == 2:
        ## [Loop over the x_avg grid.]
        for k in range(len(x_avg)-1): 
            ## [y-vals when x-vals larger than the lower bound cell.]
            y_lbnd = y_dat[x_dat>x_avg[k]]
            x_lbnd = x_dat[x_dat>x_avg[k]]
            ## [y-vals smaller than the upper bound cell.]
            y_bnd = y_lbnd[x_lbnd<x_avg[k+1]]
            ## [Set the y-val avg]
            y_avg[k, :] = np.average(y_bnd, axis=0)
        ## [Duplicate the last entry of y_avg.]
        y_avg[-1,:] = y_avg[-2,:]


    return y_avg


def Smooth_Profile_55(x_dat, y_dat): 
    """
    This function implements the 557 data smoothing algorithm as described and used by 
    Chris Wingard, but without the 7. 

    Parameters
    ----------
    x_dat: 1-D Array, [N]
        The x-coordinates of the data.
    y_dat: 1-D Array, [N]
        The y-coordinates of the data.

    Returns 
    -------
    x_55: 1-D Array, [int(N)/5/5]
        The returned x_coordinates of the smooth y data.
    y_55: 1-D Array, [int(N)/5/5]
        The returned smoothed y_data.
    """

    N = len(y_dat)

    ## [5 point median section.]
    ## -------------------------
    ## [The 5 point median array size.]
    N_5 = int(N/5)
    y_5 = np.zeros(N_5)
    x_5 = np.zeros(N_5)

    ## [Counter.]
    i = 0
    ## [Loop over profile in 5 point increments and take the 5 point median.]
    for k in range(5, N, 5): 
        ## [5 point median. Since the 5 pt array is always odd, we can assume that the 
        ##  median will be in the array.]
        y_5[i] = np.median(y_dat[k-5:k])
        ## [The correct x coordinate.]
        xi = np.argwhere(y_dat[k-5:k] == y_5[i])
        x_5[i] = x_dat[xi[0] + (k-5)]
        
        i+=1

    ## [5 point mean section.]
    ## -----------------------
    ## [The 5 point mean array size.]
    N_55 = int(N_5/5)
    y_55 = np.zeros(N_55)
    x_55 = np.zeros(N_55)

    ## [Counter.]
    i=0
    ## [Loop over the median profile in increments of 5.]
    for k in range(5, N_5, 5): 
        ## [The 5 point mean.]
        y_55[i] = np.mean(y_5[k-5:k])
        ## [x cooirdinate is the middle of the center coordinate of this 5 point mean.]
        x_55[i] = x_5[k-3]

        i+=1

    return x_55, y_55


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


def Plot_OOI_Profile(prof_index, data_set, var_name, site, assembly, instrument, method, lam=None): 
    """

    """

    ## [Read in the diffferent variables for spkir data set.]
    depth_dat = data_set.variables['depth']
    dt_dat = data_set.variables['time']
    var_dat = data_set.variables[var_name]

    dz_max = 1

    ## [Get the wavelengths for optaa]
    if var_name == 'optical_absorption': 
        assert lam != None
        wavelength_dat = data_set.variables['wavelength_a']
        depth_profs, dt_profs, wavelength_profs = Create_Profiles(depth_dat.data, dt_dat.data, wavelength_dat.data, dz_max)
        depth_profs, dt_profs, var_profs = Create_Profiles(depth_dat.data, dt_dat.data, var_dat.data, dz_max)
        abs_wavelengths = wavelength_profs[prof_index][0,:]
        lam_i = Get_Wavelength_Index(abs_wavelengths, lam)
        depth = -depth_profs[prof_index]
        dt = dt_profs[prof_index]
        var = np.squeeze(var_profs[prof_index][:,lam_i])

    else: 
        ## [Get the corresponding data in profile list form.]
        depth_profs, dt_profs, var_profs = Create_Profiles(depth_dat.data, dt_dat.data, var_dat.data, dz_max)
        ## [Make the given profs theri own variables.]
        depth = -depth_profs[prof_index]
        dt = dt_profs[prof_index]
        var = var_profs[prof_index]

    ## [Make the 1m avg grid.]
    x_avg = np.arange(depth[0], depth[-1], 1)

    ## [Plot the profile for th egiven prof_index.]
    fig, ax = plt.subplots()
    ax.plot(var, depth, '--', lw=0.5)
    ax.plot(Grid_Average_Profile(depth, var, x_avg), x_avg)
    ## [Labels.]
    ax.set_ylabel(f"Z [{depth_dat.attrs['units']}]")
    ax.set_xlabel(f"{var_name} [{var_dat.attrs['units']}]")
    ## [The date time array in seconds.]
    dt_sec = dt.astype('datetime64[s]')
    ax.set_title(f"OOI {var_name} from \n {dt_sec[0]} to {dt_sec[-1]}")
    ## [Putting some identifying text on the figure.]
    ## [10% up the vertical location]
    txt_y = ax.get_ylim()[1] + 0.5 * ax.get_ylim()[0] 
    ## [10% of the horizontal location.]
    txt_x = ax.get_xlim()[0] + 0.7 * ax.get_xlim()[1]
    ## [The change in txt location in vertical.]
    txt_dz = 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ## [Adding the txt.]
    ax.text(txt_x, txt_y, f'SITE: {site}')   
    ax.text(txt_x, txt_y+txt_dz, f'ASSEMBLY: {assembly}')   
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: {instrument}')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   
    if lam: 
        ax.text(txt_x, txt_y+4*txt_dz, f'WAVELENGTH: {abs_wavelengths[lam_i]}')   

    ax.grid()

    fig.show()

    return 


def Plot_Hoffmuller(data_set, var_name, site, assembly, instrument, method): 
    """

    """

    ## [Read in the diffferent variables for spkir data set.]
    depth_dat = data_set.variables['depth']
    dt_dat = data_set.variables['time']
    var_dat = data_set.variables[var_name]

    ## [Get the corresponding data in profile list form.]
    dz_max = 1
    depth_profs, dt_profs, var_profs = Create_Profiles(depth_dat.data, dt_dat.data, var_dat.data, dz_max)

    ## [Create the hoffmuller diagram.]
    ## [First find the largest size profile in the list of profiles.]
    ## [This largest size will set the vertical size of the mesh.]
    sizes = np.zeros(len(var_profs))
    for k,var_prof in enumerate(var_profs): 
        sizes[k] = len(var_prof)
    size = max(sizes)
        
    ## [Constructing the mesh.]
    var_mesh = np.NaN * np.zeros((int(size), len(var_profs)))
    ## [The start times of each prof for x-coords.]
    t_coord = np.zeros(len(dt_profs))

    ## [Filling the mesh with profiles, notes they are note interpolated
    ##  so this method is incorrect, but it will work for now.]
    for k, var_prof in enumerate(var_profs):
        var_len = len(var_prof)
        var_mesh[:var_len, k] = var_prof
        t_coord[k] = dt_profs[k][0] 
        if len(var_prof) == size: 
            z_coord =  depth_profs[k]

    ## [Plotting the Hovmuller disgram.]
    fig, ax=plt.subplots()
    tt, zz = np.meshgrid(t_coord, z_coord)
    im = ax.pcolormesh(tt, zz, var_mesh, cmap='nipy_spectral')
    fig.colorbar(im, ax=ax)

    fig.show()

    return 


