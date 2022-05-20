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
import xarray as xr
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


def Download_OOI_Data(ooi_savefile_head, download, site_name, assembly, method, start, stop, profiled=True):
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
    profiled: optional, Boolean
        This flag is set true if the data is to made into profiles. 

    Returns
    -------
    flort_dat: ooi data_request object
        The flort data.
    spkir_dat: ooi data_request object
        The spkir data.
    optaa_dat: ooi data_request object
        The optaa data.
    ---
    AND 
    ---

    flort_profiles: List of ooi data_request objects
        The flort data.
    spkir_profiles: List of ooi data_request objects
        The spkir data.
    optaa_profiles: List of ooi data_request objects
        The optaa data.
    """

    ## [The savefile names for each ooi data set.]
    spkir_savefile = ooi_savefile_head + '_spkir.p'
    flort_savefile = ooi_savefile_head + '_flort.p'
    optaa_savefile = ooi_savefile_head + '_optaa.p'
    ## [The pickle file names of the processed data.]
    spkir_profiles_savefile = ooi_savefile_head + '_spkir_profiles.p'
    flort_profiles_savefile = ooi_savefile_head + '_flort_profiles.p'
    optaa_profiles_savefile = ooi_savefile_head + '_optaa_profiles.p'

    ## [Download OOI data.]
    if download: 
        ## [Download the data.]
        ## [The Chla data.]
        flort_dat = Download_Data(site_name, assembly, 'FLORT', method, start, stop)
        ## [The spkir data.]
        spkir_dat = Download_Data(site_name, assembly, 'SPKIR', method, start, stop)
        ## [The OPTAA data.]
        optaa_dat = Download_Data(site_name, assembly, 'OPTAA', method, start, stop)

        ## [Make into profiles if flag True.]
        if profiled: 
            ## [Making Flort profiles, no profile processing for now.]
            flort_profiles = Create_Profiles(flort_dat, process_profile=False)
            ## [Spkir profiles.]
            spkir_profiles = Create_Profiles(spkir_dat, process_profile=False)
            ## [The optaa profiles, does include processing since the optaa data set 
            ##  has the on_seconds variable.]
            optaa_profiles = Create_Profiles(optaa_dat, process_profile=False)

            ## [Save the profiled data into pickles.]
            pickle.dump(flort_profiles, open(flort_profiles_savefile, 'wb'))
            pickle.dump(spkir_profiles, open(spkir_profiles_savefile, 'wb'))
            pickle.dump(optaa_profiles, open(optaa_profiles_savefile, 'wb'))

        ## [Save the ooi data into the pickle file.]
        pickle.dump(flort_dat, open(flort_savefile, 'wb'))
        pickle.dump(spkir_dat, open(spkir_savefile, 'wb'))
        pickle.dump(optaa_dat, open(optaa_savefile, 'wb'))



    ## [Load the data from pickle.]
    else:
        flort_dat = pickle.load(open(flort_savefile, 'rb'))
        spkir_dat = pickle.load(open(spkir_savefile, 'rb'))
        optaa_dat = pickle.load(open(optaa_savefile, 'rb'))
        if profiled: 
            flort_profiles = pickle.load(open(flort_profiles_savefile, 'rb'))
            spkir_profiles = pickle.load(open(spkir_profiles_savefile, 'rb'))
            optaa_profiles = pickle.load(open(optaa_profiles_savefile, 'rb'))

    if profiled: 
        return flort_dat, spkir_dat, optaa_dat, flort_profiles, spkir_profiles, optaa_profiles
    else:
        return flort_dat, spkir_dat, optaa_dat

def Process_Profile(profile, drop_time=60): 
    """
    This function processes a given profile. 

    Much of this function comes from line 41 of the following Jupyter notebook from Chris Wingard: 
    https://nbviewer.org/github/cwingard/ooi-data-explorations/blob/notebooks/python/examples/notebooks/optaa/process_ooinet_optaa_cspp.ipynb

    Parameters
    ----------
    profile:  OOI xarray.DataArray object
        The profiled data set. 
    drop_time: optional, Float
        The default is 60[s]. The first number of seconds of the profile to drop. 

    Returns
    -------
    profile: OOI xarray.DataArray object
        The processed profile.
    """

    ## [First remove the first 60 seconds of the profile. This is because the filter wheel spin-up and 
    ##  the lamp warmup occur during this time.]
    profile = profile.where(profile['on_seconds'] > drop_time, drop=True)

    ## [Bin the data into 25cm depth bins.]
    ## [The bins centers.]
    bin_centers = np.arange(0.125, 75.125, 0.25)
    ## [This results in a GroupBy Object with key, value pairs for each group.]
    bins = profile.groupby_bins('depth', bin_centers) 

    binned = []
    ## [Loop over each bin. Each bin is a group of a interval key and a xarray.DataArray object value, 
    ##  thus each group will be indexed with grp[1], for the group value.]
    for grp in bins: 
        ## [Taking the temporal median of each bin.]
        med = grp[1].median('time', keepdims=True, keep_attrs=True)
        ## [To get back the time coordinate is necessary to use the mean of the 
        ##  time coordinates of this group. This is because the median function does not 
        ##  return the coordinate its is taken of.  This is because the median does not 
        ##  necessarily correspond to an actual point in the xarray.DataArray object.]
        med = med.assign_coords({'time': np.atleast_1d(grp[1]['time'].mean().values)})
        ## [Set the depth coordinate to the bin midpoint.]
        med['depth'] = med['depth']*0 + grp[0].mid
        ## [Append the median data set back into a list of bins.]
        binned.append(med)

    ## [Recombine the data into a single data set.]
    profile = xr.concat(binned, 'time')  
    ## [Sort ascending in time.]
    profile = profile.sortby('time')
        
    return profile


def Create_Profiles(ooi_dat, dt_sep=120, process_profile=True): 
    """
    This organizes the x_array data set into a list of different profiles.

    Much of this function comes from line 42 of the following Jupyter notebook from Chris Wingard: 
    https://nbviewer.org/github/cwingard/ooi-data-explorations/blob/notebooks/python/examples/notebooks/optaa/process_ooinet_optaa_cspp.ipynb

    Parameters
    ----------
    ooi_data:  OOI xarray.DataArray object
        The data set. 
    dt_sep: optional, Float
        The default is 120 [s]. The minimum time between profiles in seconds.
    process_profile: optional, Boolean
        The default is True. This is used to implement the profile processing function on each profile.

    Returns
    -------
    profiles: List
        A list of profiled xarrays.
    """

    ## [This gets the index of each profile split.]
    dt = ooi_dat.where(ooi_dat['time'].diff('time') > np.timedelta64(120, 's'), drop=True).get_index('time') 

    ## [Init the list of profiled xarrays.]
    profiles = []
    ## [Process each profile by looping over the dts.]
    for k, d in enumerate(dt): 
        ## [The first profile.]
        if k==0: 
            profile = ooi_dat.where(ooi_dat['time'] < d, drop=True)
        else:
            profile = ooi_dat.where((ooi_dat['time'] >= dt[k-1]) & (ooi_dat['time'] < d), drop=True)

        ## [Append the profiled data set to the list.]
        if process_profile == True: 
            profiles.append(Process_Profile(profile))
        else: 
            profiles.append(profile)
    
    ## [Append the last profile to the list.]
    profile = ooi_dat.where(ooi_dat['time'] >= dt[-1], drop=True)
    ## [Append the profiled data set to the list.]
    if process_profile == True: 
        profiles.append(Process_Profile(profile))
    else: 
        profiles.append(profile)
                
    return profiles


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
        ## [If the length of xi is zero then it wasn't found, just use 0.]
        ## [This is a hack job fix but check with jonathan and Chris about it later.]
        if len(xi)==0: 
            x_5[i] = x_dat[0 + (k-5)]
        else:
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


