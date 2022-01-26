"""
This file is to read in the data from the plymouth laboratory ocean color data from their
OPeNDAP server, using the specific url. 
The function should take in the base url of the data set and then download data for the 
desired time period and domain
"""

## External Mods. 
import datetime
from netCDF4 import Dataset
import numpy as np
## User Mods
import cal_data 

def Days_To_Julian_Date(start_date, days_array):
    """
    Turns a the number of days since a given date into the 
    resulting date
    """
    ## Empty julian date array.
    julian_date = np.zeros_like(days_array)
    year = np.zeros(len(days_array))

    ## Loop over the array.
    for k, days in enumerate(days_array):
        date_1 = datetime.datetime.strptime(start_date, "%m/%d/%Y")
        date = date_1 + datetime.timedelta(days=int(days))
        date = date.strftime("%m/%d/%Y")
        ## using the Date_to_Julian_Date function from cal_data.py
        julian_date[k] = cal_data.Date_to_Julian_Date(date, fmt="%m/%d/%Y")
        ## Getting the year 
        year[k] = int(date.split('/')[2])

    return julian_date, year

def Get_PML_OC_Data_Set(erddap_url, year_lims, julian_date_lims, lat_lims, lon_lims):
    """
    This function returns the Rrs443, Rrs560, and chla arrays for the given 
    date_time domain and lattitude/longitude domain.
    """

    #ply_dat_url = 'https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v5.0-DAILY?lat[0:1:0],lon[0:1:0],time[0:1:0],Rrs_443[0:1:0][0:1:0][0:1:0],Rrs_560[0:1:0][0:1:0][0:1:0],chlor_a[0:1:0][0:1:0][0:1:0]' 

    ## First get the full lattitude, longitude, and time arrays to for indexing.
    ## Edit the url to include all of those arrays.
    coord_url = erddap_url.replace('lat[0:1:0]', 'lat[0:1:4319]')
    coord_url = coord_url.replace('lon[0:1:0]', 'lon[0:1:8639]')
    coord_url = coord_url.replace('time[0:1:0]', 'time[0:1:8493]')

    ## The dataset object.
    coord_ds = Dataset(coord_url)
    
    ## The Arrays for each coordinate.
    lat_coord = coord_ds.variables['lat'][:]
    lon_coord = coord_ds.variables['lon'][:]
    time_coord = coord_ds.variables['time'][:]
    
    ## Making the time coord into a date.
    ## The units of the time coord are in days since January 01, 1970. 
    julian_date, year = Days_To_Julian_Date('01/01/1970', time_coord)

    ## Now to get the indexing of the array from the 
    ## desired input limits.
    ## lattitude.
    ## The lower index
    lat_il = np.argmin(abs(lat_coord - lat_lims[0]))
    ## The upper index
    lat_iu = np.argmin(abs(lat_coord - lat_lims[1]))
    ## Longitude.
    lon_il = np.argmin(abs(lon_coord - lon_lims[0]))
    lon_iu = np.argmin(abs(lon_coord - lon_lims[1]))
    ## Year
    year_il = np.argmin(abs(year - year_lims[0]))
    year_iu = np.argmin(abs(year - year_lims[1]))
    ## upper year lim mmask
    year_maskl = year == year_lims[0]
    year_masku = year == year_lims[1]
    ## julian day
    date_il = year_il + np.argmin(abs(julian_date[year_maskl] - julian_date_lims[0]))
    date_iu = year_iu + np.argmin(abs(julian_date[year_masku] - julian_date_lims[1]))

    ## Now to remake the url with the given index bounds
    bound_url = erddap_url.replace('lat[0:1:0]', f'lat[{lat_il}:1:{lat_iu}]')
    bound_url = bound_url.replace('lon[0:1:0]', f'lon[{lon_il}:1:{lon_iu}]')
    bound_url = bound_url.replace('time[0:1:0]', f'time[{date_il}:1:{date_iu}]')
    ## Now to edit the index bounds of the Rrs443, Rrs560, and chla
    ## The indexs are of the order time[0:1:0]lat[0:1:0]lon[0:1:0]
    bound_url = bound_url.replace('Rrs_443[0:1:0][0:1:0][0:1:0]', f'Rrs_443[{date_il}:1:{date_iu}][{lat_il}:1:{lat_iu}][{lon_il}:1:{lon_iu}]')
    bound_url = bound_url.replace('Rrs_560[0:1:0][0:1:0][0:1:0]', f'Rrs_560[{date_il}:1:{date_iu}][{lat_il}:1:{lat_iu}][{lon_il}:1:{lon_iu}]')
    bound_url = bound_url.replace('chlor_a[0:1:0][0:1:0][0:1:0]', f'chlor_a[{date_il}:1:{date_iu}][{lat_il}:1:{lat_iu}][{lon_il}:1:{lon_iu}]')

    ## Get the dataset object with the bound url
    pml_ds = Dataset(bound_url)
    
    return pml_ds
    
    
def Get_Point_PML_Dataset(pml_ds, year, julian_day, lat, lon):  
    """
    This function takes the pml dataset as an argument and returns the 
    Rrs443, Rrs560, and chlor_a value for the desired spatial and temporal
    locations. 
    """

    lat_pml = pml_ds.variables['lat'][:]
    lon_pml = pml_ds.variables['lon'][:]
    time_pml = pml_ds.variables['time'][:]
    
    ## The units of the time coord are in days since January 01, 1970. 
    julian_day_pml, year_pml = Days_To_Julian_Date('01/01/1970', time_pml)
    
    ## Now to get the index of the nearest neighbour in the PLM dataset. 
    lat_i = np.argmin(abs(lat_pml - lat))
    lon_i = np.argmin(abs(lon_pml - lon))
    year_i = np.argmin(abs(year_pml - year))
    year_mask = year_pml == year
    time_i = year_i + np.argmin(abs(julian_day_pml[year_mask] - julian_day))

    ## Next to splice the desired value arrays accordingly. 
    Rrs443 = pml_ds.variables['Rrs_443'][time_i, lat_i, lon_i]
    Rrs560 = pml_ds.variables['Rrs_560'][time_i, lat_i, lon_i]
    chla = pml_ds.variables['chlor_a'][time_i, lat_i, lon_i]

    return Rrs443, Rrs560, chla





















    

    


    
