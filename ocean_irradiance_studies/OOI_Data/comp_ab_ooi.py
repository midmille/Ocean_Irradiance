"""
comp_ab_ooi.py

Author: Miles D. Miller, University of California Santa Cruz

Created: 12:18pm April 13, 2022 

Purpose: The analysis and exploration of the ooi derived optical absorption and attenuation from 
         raw data product. 
"""

## [User Modules.]
import OOI_Data_Functions as ODF
## [External Modules.]
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import io
import os 
import pandas as pd 
import re 
import pickle
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
## [OOI Modules]
from ooi_data_explorations.data_request import data_request
from ooi_data_explorations.common import m2m_request, m2m_collect, get_deployment_dates
from pyseas.data.opt_functions import opt_internal_temp, opt_external_temp
from pyseas.data.opt_functions import opt_pd_calc, opt_tempsal_corr
from cgsn_processing.process.finding_calibrations import find_calibration
from cgsn_processing.process.proc_optaa import Calibrations

def Downlaod_Chris_Depl(): 
    """
    """

    ## [Setup needed parameters for the request.]
    site = 'CE02SHSP'           # OOI Net site designator
    node = 'SP001'              # OOI Net node designator
    sensor = '04-OPTAAJ000'     # OOI Net sensor designator
    method = 'recovered_cspp'   # OOI Net data delivery method
    stream = 'optaa_dj_cspp_instrument_recovered'  # OOI Net stream name
    start, stop = get_deployment_dates(site, node, sensor, 19)
    tag = '.*deployment0019.*OPTAA.*\\.nc$'  # limit request to OPTAA NetCDF files from Deployment 19

    # request the data from the OOI M2M API ...
    r = m2m_request(site, node, sensor, method, stream, start, stop)
    data = m2m_collect(r, tag)

    return data

def get_dev_coeffs(data): 
    """
    """
    # load the instrument calibration data
    coeff_file = 'optaa.cal_coeffs.json'
    dev = Calibrations(coeff_file)  # initialize calibration class
    # check for the source of calibration coeffs and load accordingly
    if os.path.isfile(coeff_file):
        # we always want to use this file if it already exists
        dev.load_coeffs()
    else:
        # load from the CI hosted CSV files
        csv_url = find_calibration('OPTAA', data.attrs['SerialNumber'][4:], 
                                        data['time'][0].values.astype(float) / 10**9)
        if csv_url:
            tca_url = re.sub('.csv', '__CC_taarray.ext', csv_url)
            tcc_url = re.sub('.csv', '__CC_tcarray.ext', csv_url)
            dev.read_devurls(csv_url, tca_url, tcc_url)
            dev.save_coeffs()
    # check the device file coefficients against the data file contents
    if dev.coeffs['serial_number'] != int(data.attrs['SerialNumber'][4:]):
        raise Exception('Serial Number mismatch between ac-s data and the device file.')
    print('num_wavelengths', dev.coeffs['num_wavelengths'], data['num_wavelengths'][0])
#    if dev.coeffs['num_wavelengths'] != data['num_wavelengths'][0]):
#        raise Exception('Number of wavelengths mismatch between ac-s data and the device file.')
                                                                                                     
    return dev


def apply_dev(optaa, coeffs):
    """
    Processes the raw data contained in the optaa dictionary and applies the 
    factory calibration coefficents contained in the coeffs dictionary to
    convert the data into initial science units.
 
    :param optaa: xarray dataset with the raw absorption and beam attenuation
    measurements.
    :param coeffs: Factory calibration coefficients in a dictionary structure
 
    :return optaa: xarray dataset with the raw absorption and beam attenuation
    measurements converted into particulate and beam attenuation values
    with the factory pure water calibration values subtracted.
    """
    # convert internal and external temperature sensors
    optaa['internal_temp'] = opt_internal_temp(optaa['internal_temp_raw'])
    optaa['external_temp'] = opt_external_temp(optaa['external_temp_raw'])
 
    # setup inputs
    a_ref = optaa['a_reference_counts']
    a_sig = optaa['a_signal_counts']
    c_ref = optaa['c_reference_counts']
    c_sig = optaa['c_signal_counts']
    npackets = a_ref.shape[0]
 
    # initialize the output arrays
    apd = a_ref * np.nan
    cpd = c_ref * np.nan
     
    # calculate the L1 OPTAA data products (uncorrected beam attenuation and absorbance) for particulate
    # and dissolved organic matter with pure water removed.
    for ii in range(npackets):
        # calculate the uncorrected optical absorption coefficient [m^-1]
        apd[ii, :], _ = opt_pd_calc(a_ref[ii, :], a_sig[ii, :], coeffs['a_offsets'],
                                    optaa['internal_temp'].values[ii], coeffs['temp_bins'],
                                    coeffs['ta_array'])
        # calculate the uncorrected optical attenuation coefficient [m^-1]
        cpd[ii, :], _ = opt_pd_calc(c_ref[ii, :], c_sig[ii, :], coeffs['c_offsets'],
                                    optaa['internal_temp'].values[ii], coeffs['temp_bins'],
                                    coeffs['tc_array'])
 
    # save the results back to the data set
    optaa['apd'] = apd
    optaa['cpd'] = cpd
 
    # return the optaa dictionary with the factory calibrations applied
    return optaa
 
 
def apply_tscorr(optaa, coeffs, temp=None, salinity=None):
    """
    Corrects the absorption and beam attenuation data for the absorption
    of seawater as a function of the seawater temperature and salinity (the
    calibration blanking offsets are determined using pure water.)
     
    If inputs temp or salinity are not supplied as calling arguments, then the 
    following default values are used.
     
    temp: temperature values recorded by the ac-s's external thermistor.
    salinity: 33.0 psu
 
    Otherwise, each of the arguments for temp and salinity should be either a 
    scalar, or a 1D array or a row or column vector with the same number of time
    points as 'a' and 'c'.
 
    :param optaa: xarray dataset with the temperature and salinity corrected
    absorbance data array that will be corrected for the effects of
    scattering.
    :param coeffs: Factory calibration coefficients in a dictionary structure
    :param temp: In-situ seawater temperature, ideally from a co-located CTD
    :param salinity: In-situ seawater salinity, ideally from a co-located CTD
 
    :return optaa: xarray dataset with the temperature and salinity corrected
    absorbance and attenuation data arrays added.
    """
    # setup the temperature and salinity arrays
    if temp is None:
        temp = optaa['external_temp'].values
    else: 
        if np.array(temp).size == 1:
            temp = np.ones_like(optaa['external_temp']) * temp
        else:
            temp = np.array(temp)
     
    if temp.size != optaa['time'].size:
        raise Exception("Mismatch: temperature array != number of OPTAA measurements")
 
    if salinity is None:
        salinity = np.ones_like(optaa['external_temp']) * 33.0
    else:
        if np.array(salinity).size == 1:
            salinity = np.ones_like(optaa['external_temp']) * salinity
        else:
            salinity = np.array(salinity)
 
    if salinity.size != optaa['time'].size:
        raise Exception("Mismatch: salinity array != number of OPTAA measurements")
 
    # setup and size the inputs
    apd = optaa['apd']
    cpd = optaa['cpd']
    npackets = apd.shape[0]
 
    # initialize the output arrays
    apd_ts = apd * np.nan
    cpd_ts = cpd * np.nan
 
    # apply the temperature and salinity corrections
    for ii in range(npackets):
        wvlngth = ~np.isnan(apd[ii, :].values)
        apd_ts[ii, wvlngth] = opt_tempsal_corr('a', apd[ii, wvlngth], coeffs['a_wavelengths'][wvlngth],
                                                coeffs['temp_calibration'], temp[ii], salinity[ii])
        cpd_ts[ii, wvlngth] = opt_tempsal_corr('c', cpd[ii, wvlngth], coeffs['c_wavelengths'][wvlngth],
                                                coeffs['temp_calibration'], temp[ii], salinity[ii])
     
    # save the results
    optaa['apd_ts'] = apd_ts
    optaa['cpd_ts'] = cpd_ts
    return optaa
 
 
def apply_scatcorr(optaa, coeffs):
    """
    Correct the absorbance data for scattering using Method 1, with the
    wavelength closest to 715 nm used as the reference wavelength for the
    scattering correction.
 
    :param optaa: xarray dataset with the temperature and salinity corrected
    absorbance data array that will be corrected for the effects of
    scattering.
    :param coeffs: Factory calibration coefficients in a dictionary structure
 
    :return optaa: xarray dataset with the method 1 scatter corrected
    absorbance data array added.
    """
    # find the closest wavelength to 715 nm
    reference_wavelength = 715.0
    idx = np.argmin(np.abs(coeffs['a_wavelengths'] - reference_wavelength))
 
    # use that wavelength as our scatter correction wavelength
    apd_ts = optaa['apd_ts']
    apd_ts_s = apd_ts - apd_ts[:, idx]
 
    # save the results
    optaa['apd_ts_s'] = apd_ts_s
    return optaa
 
 
def estimate_chl_poc(optaa, coeffs):
    """
    Derive estimates of Chlorophyll-a and particulate organic carbon (POC)
    concentrations from the temperature, salinity and scatter corrected
    absorption and beam attenuation data.
 
    :param optaa: xarray dataset with the scatter corrected absorbance data.
    :param coeffs: Factory calibration coefficients in a dictionary structure
 
    :return optaa: xarray dataset with the estimates for chlorophyll and POC
    concentrations added.
    """
    # use the standard chlorophyll line height estimation with an extinction coefficient of 0.020.
    m676 = np.argmin(np.abs(coeffs['a_wavelengths'] - 676.0))
    m650 = np.argmin(np.abs(coeffs['a_wavelengths'] - 650.0))
    m715 = np.argmin(np.abs(coeffs['a_wavelengths'] - 715.0))
    apg = optaa['apd_ts_s']
    aphi = apg[:, m676] - 39/65 * apg[:, m650] - 26/65 * apg[:, m715]
    optaa['estimated_chlorophyll'] = aphi / 0.020
 
    # estimate the POC concentration from the attenuation at 660 nm
    m660 = np.argmin(np.abs(coeffs['c_wavelengths'] - 660.0))
    cpg = optaa['cpd_ts']
    optaa['estimated_poc'] = cpg[:, m660] * 380
 
    print('here')
    return optaa


def Process_Raw_Optaa(optaa_dat): 
    """
    This function implements the functions provided by Chris Wingard to show the processing 
    of raw optaa data to get final absorption adn attenuation products.
    Much of this function comes from the following Jupyter notebook from Chris Wingard: 
    https://nbviewer.org/github/cwingard/ooi-data-explorations/blob/notebooks/python/examples/notebooks/optaa/process_ooinet_optaa_cspp.ipynb
 
    Parameters
    ----------

    Returns
    -------
    """

    ## [Get the dev coefficients.]
    dev = get_dev_coeffs(optaa_dat)

    ## [Applying the functions to the data set.]
    optaa_dat = apply_dev(optaa_dat, dev.coeffs)
    optaa_dat = apply_tscorr(optaa_dat, dev.coeffs, optaa_dat.temperature, optaa_dat.salinity)
    optaa_dat = apply_scatcorr(optaa_dat, dev.coeffs)

    return optaa_dat
    
def Plot_Comp_Ab(optaa_profs, ip): 
    """
    Plots the absorption and scattering as functions of wavelength.
     
    """
    
    fig = plt.figure()
    ## [Grid spec object for subplots()] 
    gs = mpl.gridspec.GridSpec(3, 2)
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[1,1])
    ## [The axs for the absorption and scattering spectral plots.]
    axs = [ax0, ax1, ax2, ax3]
    ## [The ax for informative text.]
    ax_txt = fig.add_subplot(gs[2,:])

    
    optaa_prof = optaa_profs[ip]
#    optaa_prof = optaa_prof.dropna('wavelength')

    for k in range(optaa_prof['depth'].shape[0]): 
        ## [The different processed variables.]
        a = optaa_prof['optical_absorption'].data[k,:]
        a_mask = ~np.isnan(a)
        a = a[a_mask]
        ## !!!!!!!!IMPORTANT THIS WAVELENGTH SPLICING CHANGE IS ONLY FOR m2m download!!!!!!!
        a_lam = optaa_prof['wavelength_a'].data[k,:][a_mask]
        c = optaa_prof['beam_attenuation'].data[k,:]
        c_mask = ~np.isnan(c)
        c = c[c_mask]
        ## !!!!!!!!IMPORTANT THIS WAVELENGTH SPLICING CHANGE IS ONLY FOR m2m download!!!!!!!
        c_lam = optaa_prof['wavelength_c'].data[k,:][c_mask]
        ## [Interpolate the attenuation to the absorption grid.]
        c = np.interp(a_lam, c_lam, c) 
        ## [Scattering is simply the attenuation minus the absorption.]
        b = c - a
 
        ## [The different raw variables.]
        a_raw = optaa_prof['apd_ts_s'].data[k,:][a_mask]
        c_raw = optaa_prof['cpd_ts'].data[k,:][c_mask]
        ## [Interpolate the attenuation to the absorption grid.]
        c_raw = np.interp(a_lam, c_lam, c_raw) 
        ## [Scattering is simply the attenuation minus the absorption.]
        b_raw = c_raw - a_raw
 
        ## [Plotting the different absorption and scattering.]
        axs[0].plot(a_lam, a, 'b', linewidth = 0.7 )
        axs[1].plot(a_lam, b, 'b', linewidth = 0.7 )
        axs[2].plot(a_lam, a_raw, 'b', linewidth = 0.7 )
        axs[3].plot(a_lam, b_raw, 'b', linewidth = 0.7 )

    ## [Loop over the axes for the labels.]
    for iax, ax in enumerate(axs): 

        ## [The labels and other information.]
        ## [Labels.]
        if (iax==2) or (iax==3):
            ax.set_xlabel(f"Wavelength [{optaa_prof.variables['wavelength_a'].attrs['units']}]")
        if (iax==0) or (iax==2):
            ax.set_ylabel(f"Absorption [{optaa_prof.variables['optical_absorption'].attrs['units']}]")
        if (iax==1) or (iax==3):
            ax.set_ylabel(f"Scattering [{optaa_prof.variables['optical_absorption'].attrs['units']}]")
        ## [The titles.]
        if iax == 0: 
            ax.set_title('optical_absorption')
        if iax == 1: 
            ax.set_title('scattering = beam_attenuation - optical_absorption')
        if iax == 2: 
            ax.set_title('apd_ts_s')
        if iax == 3: 
            ax.set_title('scattering = cpd_ts - apd_ts_s')

        ax.grid()

    ## [The ax for the info.]
#    dt_sec = dt_profs[prof_index].data.astype('datetime64[s]')
    ## [Putting some identifying text on the figure.]
    ## [10% up the vertical location]
    txt_y = -0.1
    ## [10% of the horizontal location.]
    txt_x = 0.0
    ## [The change in txt location in vertical.]
    txt_dz = 0.15
    ## [Removing the ax frame.]
    ax_txt.axis('off')
    ## [Adding the txt.]
    ax_txt.text(txt_x, txt_y, f"Site: {optaa_prof.attrs['subsite']}")   
    ax_txt.text(txt_x, txt_y+txt_dz, f"Time Start: {optaa_prof.attrs['time_coverage_start']}")   
    ax_txt.text(txt_x, txt_y+2*txt_dz, f"Profile Number: {prof_index}")   
    ax_txt.text(txt_x, txt_y+3*txt_dz, f"Stream: {optaa_prof.attrs['stream']}")   
    ax_txt.text(txt_x, txt_y+4*txt_dz, f"Sensor: {optaa_prof.attrs['sensor']}")   
    ax_txt.text(txt_x, txt_y+5*txt_dz, "Details: A Figure Showing the spectraly dependent absorbtion and scattering for all depths in a single profile.")

    fig.show()
    
    return 
     
    
if __name__ == '__main__': 
 
    ## [The data request parameters.]
    site_name = "CE02SHSP"
    assembly = "profiler"
    method = "recovered_cspp"
    start = "2021-04-01"
    stop = "2021-04-30"
 
    ## [File name parameter.]
    ooi_savefile_head = "ooi_data/ooi_dat"
 
    ## [Data load flag, fale means load data from pickle.]
    download = False
    ## [Process the data and save processed data if True. Load if False.]
    process_raw = True
    ## [Use the deployement provided by Chris Wingard in his Jupyter notebook]
    use_chris_depl = True
 
    ## [Functional parameters.]
    prof_index=1
    wavelengths = np.arange(425, 725, 25)
    #    phy_species = ['HLPro', 'LLPro', 'Cocco', 'Diat', 'Syn', 'Tricho', 'Lgeuk'] 
    phy_species = ['HLPro', 'Cocco', 'Diat', 'Syn'] 
    cdom_reflam = 412.0
 
    ## [Download or load the data sets either from Chris' deployement or mine.]
    if use_chris_depl: 
        chris_savefile = "ooi_data/chris_ooi_depl.p"
        if os.path.exists(chris_savefile): 
            optaa_dat = pickle.load(open(chris_savefile, 'rb'))
        else: 
            optaa_dat = Downlaod_Chris_Depl()
            pickle.dump(optaa_dat, open(chris_savefile, 'wb'))
    else:
        ooi_data = ODF.Download_OOI_Data(ooi_savefile_head, 
                                        download, 
                                        site_name, 
                                        assembly, 
                                        method, 
                                        start, 
                                        stop)
 
        flort_dat, spkir_dat, optaa_dat, flort_profs, spkir_profs, optaa_profs = ooi_data
                                                         
        ## [Index the profiles for the given index.]
        flort_prof = flort_profs[prof_index]
        optaa_prof = optaa_profs[prof_index]
 
    ## [processs data flag.]
    if process_raw:
        ## [Processing the optaa raw part.] 
        optaa_dat = Process_Raw_Optaa(optaa_dat)
    
        ## [Remake profiles for raw optaa_dat.]
        optaa_profs = ODF.Create_Profiles(optaa_dat, process_profile=True) 
     
        ## [Save or load:]
        pickle.dump(optaa_profs, open(ooi_savefile_head+'_rawoptaaprofs.p', 'wb'))

    else: 
        optaa_profs = pickle.load(open(ooi_savefile_head+'_rawoptaaprofs.p', 'rb'))

    ## [plot the raw and final absorption products.]
    Plot_Comp_Ab(optaa_profs, prof_index)
