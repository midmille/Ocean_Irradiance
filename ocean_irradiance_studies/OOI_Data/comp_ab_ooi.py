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
import csv
import io
import os 
import pandas as pd 
import re 
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
        if dev.coeffs['num_wavelengths'] != data['num_wavelengths'][0]:
            raise Exception('Number of wavelengths mismatch between ac-s data and the device file.')
                                                                                                       
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
        wvlngth = ~np.isnan(apd[0, :].values)
     
        # initialize the output arrays
        apd_ts = apd * np.nan
        cpd_ts = cpd * np.nan
     
        # apply the temperature and salinity corrections
        for ii in range(npackets):
            apd_ts[ii, wvlngth] = opt_tempsal_corr('a', apd[ii, wvlngth], coeffs['a_wavelengths'],
                                                    coeffs['temp_calibration'], temp[ii], salinity[ii])
            cpd_ts[ii, wvlngth] = opt_tempsal_corr('c', cpd[ii, wvlngth], coeffs['c_wavelengths'],
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
    
    ## [Get the dev coefficients.]
    dev = get_dev_coeffs(optaa_dat)

    ## [Applying the functions to the data set.]
    optaa_dat = apply_dev(optaa_dat, dev.coeffs)
    optaa_dat = apply_tscorr(optaa_dat, dev.coeffs, optaa_dat.temperature, optaa_dat.salinity)
    optaa_dat = apply_scatcorr(optaa_dat, dev.coeffs)

    return optaa_dat
    
def Comp_Ab(): 
    """
     
    """
     
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
 
    ## [Functional parameters.]
    prof_index=1
    wavelengths = np.arange(425, 725, 25)
    #    phy_species = ['HLPro', 'LLPro', 'Cocco', 'Diat', 'Syn', 'Tricho', 'Lgeuk'] 
    phy_species = ['HLPro', 'Cocco', 'Diat', 'Syn'] 
    cdom_reflam = 412.0
 
    ## [Download or load the data sets.]
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
 
    ## [Processing the optaa raw part.] 
    optaa_dat = Process_Raw_Optaa(optaa_dat)
    
