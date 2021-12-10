"""
viirs_calcofi_val.py
 
Created: 10 December 2021

Author: Miles Miller


This file is to create a mathcup of the calcofi in situ data to the viirs satellite data. 

"""

## External Mods
from netCDF4 import Dataset
import numpy as np

def find_viirs_Rrs(julian_day, year, lat, lon): 
    """
    Find the nearest neighbour point of chla in the viirs domain. 
    

    Parameters
    ----------
    julian_day: Integer
        The day of the year out of 365.
    year: Integer
        The given year of the data point.
    lat: Float 
        The lattitude of the desired point.
    lon: Float 
        The longitutde of the desired point. 

    Returns
    -------
 
    """
    
    ## F0 is the conversion factor from nLw to Rrs
    ## See john Wilkin email.
    F0443 = 190.2609
    F0551 = 184.2347

    ## The normalized water leaving radiance
    nlw443_df = Dataset(f'V{year}{julian_day}_D1_WW00_nLw443','r')
    nlw551_df = Dataset(f'V{year}{julian_day}_D1_WW00_nLw551','r')
  
    ## Getting the nearest lon, lat indexes
    ## 443
    lat_443 = nlw443_df.variables['lat'][:]
    lon_443 = nlw443_df.variables['lon'][:]
    lat_i_443 = np.where( abs(lat_443 - lat) == np.minimum(abs(lat_443-lat)))
    lon_i_443 = np.where( abs(lon_443 - lon) == np.minimum(abs(lon_443-lon)))
    ## 551
    lat_551 = nlw551_df.variables['lat'][:]
    lon_551 = nlw551_df.variables['lon'][:]
    lat_i_551 = np.where( abs(lat_551 - lat) == np.minimum(abs(lat_551-lat)))
    lon_i_551 = np.where( abs(lon_551 - lon) == np.minimum(abs(lon_551-lon)))

    ## nLw 
    nlw443_arr = nlw443_df.variables['nLw_443'][:]
    nlw551_arr = nlw551_df.variables['nLw_551'][:]

    nlw443 = nlw443_arr[0,0,lat_i_443, lon_i_443]
    nlw551 = nlw551_arr[0,0,lat_i_551, lon_i_551]

    ## Making into Rrs.
    Rrs443 = nlw443 / F0443    
    Rrs551 = nlw551 / F0551

    return Rrs443, Rrs551




