"""
viirs_calcofi_val.py
 
Created: 10 December 2021

Author: Miles Miller


This file is to create a mathcup of the calcofi in situ data to the viirs satellite data. 

"""
## User Mods
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR

## External Mods
from netCDF4 import Dataset
import numpy as np
import os

def Viirs_Rrs_Chla(julian_day, year, lat, lon): 
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

    ## Padding day with zeros.
    if julian_day < 10:
        julian_day = f'00{julian_day}'
    elif julian_day < 100: 
        julian_day = f'0{julian_day}'
    else: 
        julian_day = f'{julian_day}'

    path_443 = f'VIIRS_Downloads/V{year}{julian_day}_D1_WW00_nLw443.nc'
    path_551 = f'VIIRS_Downloads/V{year}{julian_day}_D1_WW00_nLw551.nc'
    ## The normalized water leaving radiance
    ## checking to ensure the download file exists
    if os.path.exists(path_443) == False or os.path.exists(path_551) == False:
        print('Path does not exist')
        print(path_443)
        print(path_551)
        return np.NAN, np.NAN, np.NAN, np.NAN, np.NAN
    
    nlw443_df = Dataset(path_443,'r')
    nlw551_df = Dataset(path_551,'r')
  
    ## nLw 
    nlw443_arr = nlw443_df.variables['nLw_443'][0,0,:,:]
    nlw551_arr = nlw551_df.variables['nLw_551'][0,0,:,:]
    
    ## Getting the nearest lon, lat indexes
    ## 443
    lat_443 = nlw443_df.variables['lat'][:]
    lon_443 = nlw443_df.variables['lon'][:]
    lat_i_443 = np.argmin( abs(lat_443 - lat) )
    print(lat_i_443)
    lon_i_443 = np.argmin( abs(lon_443 - lon))

    ## 551
#    lat_551 = nlw551_df.variables['lat'][:]
#    lon_551 = nlw551_df.variables['lon'][:]
#    lat_i_551 = np.argmin( abs(lat_551 - lat) )
#    lon_i_551 = np.argmin( abs(lon_551 - lon) )

    near_lat = lat_443[lat_i_443]
    near_lon = lon_443[lon_i_443]

    nlw443 = nlw443_arr[lat_i_443, lon_i_443]
    nlw551 = nlw551_arr[lat_i_443, lon_i_443]

    ## Making into Rrs.
    Rrs443 = nlw443 / F0443    
    Rrs551 = nlw551 / F0551

    ## calculating chla from rrs
    chla = OIR.OCx_alg(Rrs443, Rrs551)
    

    return Rrs443, Rrs551, chla, near_lat, near_lon




