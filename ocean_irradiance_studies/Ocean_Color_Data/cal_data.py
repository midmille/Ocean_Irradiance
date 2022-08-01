"""
Created: November 14 2021 

Author: Miles Miller


This file is to automate the data analysis of the NASA ocean color NOMAD data set to 
better analyze my three stream irradiance solution. 

Implementation
--------------

1. First read in the chlorophyll validation data set from the NOMAD link on the ocean color website.
   Also available in email with Jonathan.

2. The validation csv file contains the surface chlorophyll data from satllite using rrs --> chla ocean color algorithim 
   as well as in situ data taken from cruises. The validation data set just has the surface insitu but it does have 
   the name of the cruise the data was taken. I need the entire chlorophyll profile for my irradiance algorithim and that
   is where the other data sets come in. 

3. For the given validatation row there is a corresponding cruise ID. I have data from the cal cofi cruises from Jonathan. 
   Thus for the given cruise ID there is a correspnding 'cast' from the cal_cast csv file. From this file I can find the cast count
   corresponding to the cruise ID and the desired date/time of the validation data row.  

4. Then I can use that cast count and gather the chlorophyll profile corresponding to that cast count from the cal bottle csv file.

5. Then using that chlorophyll profile I can implement my irradiance algorithim and compare the resulting rrs and surface chlrophyll value.
 
"""

## External Mods
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import datetime
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import geopy.distance
## User made mod
import ocean_irradiance_module.Ocean_Irradiance as OI
#import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
#import ocean_irradiance_baird.ocean_irradiance_baird as OIB
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
import ocean_irradiance_module.Wavelength_To_RGB as W2RGB
#import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC
import viirs_calcofi_val
import cci_oc


def Get_Data(oc_chla_val_file, cal_cast_file, cal_bot_file): 
    """
    This returns a pandas data frame of the csv file from the NASA ocean color NOMAD validation data set

    """
    ## Header is row 26 and the skipped rows are units. 
    chla_val_dat = pd.read_csv(oc_chla_val_file, delimiter=',', header=26, skiprows=[27,28])   

    cal_cast_dat = pd.read_csv(cal_cast_file, delimiter=',', header=0, low_memory = False)

    cal_bot_dat =  pd.read_csv(cal_bot_file, delimiter=',', header=0, encoding='unicode_escape', low_memory=False)

    return chla_val_dat, cal_cast_dat, cal_bot_dat


def Get_Cst_Cnt(chla_val_dat, cal_cast_dat, id): 
    """
    This gets the cast count for the desired field id from a cal211 cruise. 
    """  
    
    ## boolean indexing to get df row that corresponds to id. 
    row_bool = chla_val_dat['/fields=id'] == id

    ## Then retrieve date and time from row. 
    date_time = chla_val_dat.loc[row_bool,'date_time'].item()
    date, time = date_time.split()
    
    ## rewriting the date with '/' not '-'.
    date_split = date.split('-')
    date = f'{date_split[1]}/{date_split[2]}/{date_split[0]}'
    ## Then to retireve the cast count for this given date time cal cruise.
    cast_row_date = cal_cast_dat[cal_cast_dat['Date'] == date]
    ## using the split function because the time column has a date part and a time part 
    ## and I only want the time.
    cast_row_bool = cast_row_date['Time'].str.contains(time, regex=False)
    ## Getting the cast count.
    cast_count_row = cast_row_date.loc[cast_row_bool,'Cst_Cnt']
    if cast_count_row.empty:
        cast_count = None
    else:
        cast_count = cast_count_row.item()

    return cast_count


def Get_Cast_Depth_Chla(cal_bot_dat, cast_count, incld_salt=False): 
    """
    This function returns the chla profile and z-coordinates for the given cast_count
    """

    ## This is the different bottles (points in profile) of the desired cast.
    cal_bot_set = cal_bot_dat[cal_bot_dat['Cst_Cnt'] == cast_count]
    ## The chl a profile.
    chla = cal_bot_set['ChlorA'].to_numpy()
    ## The depth profile.
    z = cal_bot_set['Depthm'].to_numpy()

    ## Formatting the profile is ROMS coordinates, bot at index 0.
    chla = np.flip(chla)
    ## making z-coordinates negative also. 
    z = -np.flip(z)

    ## Getting rid of nans
    mask = ~np.isnan(chla)
    z = z[mask]
    chla = chla[mask]

    if incld_salt: 
        salt = cal_bot_set['Salnty'].to_numpy()
        salt = np.flip(salt)
        salt = salt[mask]

        return z, chla, salt

    else:
        return z, chla

 
def Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=None): 
    """
    This function loops over the cruise data and calculates the irradiances and chla at the surface for each. 
    """ 



    ## using my params
    PI = Param_Init()

    Eu_surf = np.zeros(len(chla_val_cal_dat['cruise']))
    chla_dat = np.zeros(len(chla_val_cal_dat['cruise']))
    Eu_surf_dict = {}
    ## Looping over ocean color wavelengths.
    for lam in [443, 551]:
        ## Then loop over cal cruises. 
        for k, id in enumerate(chla_val_cal_dat['/fields=id']): 
            ## Getting the cast count. 
            cst_cnt = Get_Cst_Cnt(chla_val_cal_dat, cal_cast_dat, id)
            if cst_cnt == None: 
                Eu_surf[k] = np.nan
                continue 
            ## Getting the chla and depth profile.
            z, chla, salt = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt, incld_salt=True)     
            chla_dat[k] = chla[-1]
            ## Calculating the irradiance
            phy = OI.Phy(z, chla, ESD(phy_type), abscat(lam, phy_type, C2chla=C2chla)[0], abscat(lam, phy_type, C2chla=C2chla)[1])
            #cdom = OI.CDOM(z, salt, lam)

            irr_out = OI.ocean_irradiance(PI, 
                                          z[0], 
                                          abscat(lam, 'water'), 
                                          phy=phy, 
                                          N=1000)
            Eu_surf[k] =  irr_out[2][-1]
        Eu_surf_dict[lam] = np.copy(Eu_surf)
            
         
    ## Calculating the rrs for each wavelength and the chla 
    rrs_443 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf_dict[443])
    rrs_551 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf_dict[551])

    ## Chla 
    irr_chla = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)

    return rrs_443, rrs_551, irr_chla, chla_dat        




def Run_and_Plot_Comparison(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, species, C2chla_vals, method='shoot_up', theta_air = None):
    """
    Plots a scatter plot comparing the values
    """ 

    ## The reference data. 
    chla_insitu = chla_val_cal_dat['insitu_chlor_a']
    chla_sat = chla_val_cal_dat['aqua_chlor_a']
    rrs_443_dat = chla_val_cal_dat['aqua_rrs443'] 
    rrs_555_dat = chla_val_cal_dat['aqua_rrs555'] 
    rrs_ratio_dat =  rrs_443_dat / rrs_555_dat 

    ## This fig is for chla and rrs ratio
    fig, axes = plt.subplots(3, 2)
    ## This fig is for rrs_443, rrs_551
    fig2, axes2 = plt.subplots(3,2)
    
    ## Plot C2chla as LC2chla ratio with diff species. 
    LC2chla = 100
    for k, phy_type in enumerate(species):
        if method == 'shoot_up':
            rrs_443, rrs_551,  irr_chla, chla_dat = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=LC2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        PC.Plot_Comparison(axes[0,0], chla_insitu, irr_chla, 'Chla LC2chla Varied Species', phy_type, None, 'Model') 
        PC.Plot_Comparison(axes[0,1], rrs_ratio_dat, rrs_ratio, 'Rrs Ratio  LC2chla Varied Species', phy_type, None, None) 
        PC.Plot_Comparison(axes2[0,0], rrs_443_dat, rrs_443, 'Rrs 443 LC2chla Varied Species', phy_type, None, 'Model') 
        PC.Plot_Comparison(axes2[0,1], rrs_555_dat, rrs_551, 'Rrs 551 LC2chla Varied Species', phy_type, None, None) 
    ## plotting satelite 
    PC.Plot_Comparison(axes[0,0], chla_insitu, chla_sat, 'Chla LC2chla Varied Species', 'Sattelite', None, 'Model') 
    axes[0,0].legend( title='species')
    axes[0,0].grid()
    axes[0,1].grid()
    axes2[0,0].legend( title='species')
    axes2[0,0].grid()
    axes2[0,1].grid()

    ## Plot C2chla as SC2chla ratio with diff species. 
    SC2chla = 50
    for k, phy_type in enumerate(species):
        if method == 'shoot_up':
            rrs_443, rrs_551,  irr_chla, chla_dat = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=SC2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        PC.Plot_Comparison(axes[1,0], chla_insitu, irr_chla, 'Chla SC2chla Varied Species', phy_type, None, 'Model') 
        PC.Plot_Comparison(axes[1,1], rrs_ratio_dat, rrs_ratio, 'Rrs Ratio SC2chla Varied Species', phy_type, None, None) 
        PC.Plot_Comparison(axes2[1,0], rrs_443_dat, rrs_443, 'Rrs 443 SC2chla Varied Species', phy_type, None, 'Model') 
        PC.Plot_Comparison(axes2[1,1], rrs_555_dat, rrs_551, 'Rrs 551 SC2chla Varied Species', phy_type, None, None) 
    PC.Plot_Comparison(axes[1,0], chla_insitu, chla_sat, 'Chla SC2chla Varied Species', 'Sattelite', None, 'Model') 
    axes[1,0].legend(title='species')
    axes[1,0].grid()
    axes[1,1].grid()
    axes2[1,0].legend( title='species')
    axes2[1,0].grid()
    axes2[1,1].grid()

    ## Plot generic species with different C2chla. 
    phy_type = 'Generic'
    for k, C2chla in enumerate(C2chla_vals):
        if method == 'shoot_up':
            rrs_443, rrs_551,  irr_chla, chla_dat = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=C2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        PC.Plot_Comparison(axes[2,0], chla_insitu, irr_chla, 'Chla Generic Species Varied C2chla', C2chla, 'In Situ', 'Model') 
        PC.Plot_Comparison(axes[2,1], rrs_ratio_dat, rrs_ratio, 'Rrs Ratio Generic Species Varied C2chla', C2chla, 'In Situ', None) 
        PC.Plot_Comparison(axes2[2,0], rrs_443_dat, rrs_443, 'Rrs 443 Generic Species Varied C2chla', C2chla, 'In Situ', 'Model') 
        PC.Plot_Comparison(axes2[2,1], rrs_555_dat, rrs_551, 'Rrs 551 Generic Species Varied C2chla', C2chla, 'In Situ', None) 
    PC.Plot_Comparison(axes[2,0], chla_insitu, chla_sat, 'Chla Generic Species Varied C2chla', 'Sattelite', 'In Situ', 'Model') 
    axes[2,0].legend(title='C2chla')
    axes[2,0].grid()
    axes[2,1].grid()
    axes2[2,0].legend( title='C2chla')
    axes2[2,0].grid()
    axes2[2,1].grid()

    #plt.tight_layout()
    fig.show()
    fig2.show()
  
    ## plotting a comparison of chla dat to chla insitu 
    fig, ax = plt.subplots()
    PC.Plot_Comparison(ax, chla_insitu, chla_dat, 'Chla NOMAD In Situ Compared to Cast In Situ', None, 'NOMAD', 'Cast')
    ax.grid()
    fig.show()


    return 
 

def Run_Irr_Comp_Insitu(PI, save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, phy_type, plot=False):
    """
    This calculates the irradiance chla value for many different casts within a given time line. 
   
    Parameters
    ----------
    save_dir: String 
        The name of the directory that the irradiance field pickle file will be saved in. 
    save_file: String
        The name of the pickle file.
    wavelengths: List
        A list of the wavelengths for which irradiance will be calculated for. 
    N: Int 
        The number of vertical levels for the irradiance grid. 
    year_min: Int
        The start year for casts to be calculated. the limit is 1949. 
    cal_casts_dat: Pandas Data Frame
        The data frame containning all the casts.
    cal_bot_dat: Pandas Data Frame
        The data frame containing all the bottles corresponding to the casts.

    Returns
    -------
    irr_field: Dict
        The irradiance field dictionary where the wavelengths are keys and the field is the casts taken in more recent years
        than the year_min variable.
    """
    
 
    ## Bounded data set in the desired timeline.
    cal_cast_dat_bnd = cal_cast_dat[cal_cast_dat['Year'] > year_min]
    ## The number of casts to be calculated. 
    N_cst = len(cal_cast_dat_bnd)

    ## The ditionary of irr fields. 
    irr_field = {}

    ## Checking if the save file exists.
    ## If it doesn't then do the calculation.
    ## The location of the save file.
    save_path = f'{save_dir}/{save_file}_irr_{phy_type}.p'
    if os.path.exists(save_path) == False:
    
        for lam in wavelengths: 
            ## The array in which to store all the irradiance solutions into. 
            irr_arr = np.zeros((N,N_cst, 4))
         
            ## Now loop over the casts and calculate the irradiance chla each time.  
            for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
                print(f'{k}/{N_cst}')
                ## Getting the chla and depth profile.
                z, chla, salt = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt, incld_salt=True)     

                ## checking for zero sized return
                if len(z) == 0: 
                    z = np.zeros(N) * np.nan
                    chla = np.zeros(N) * np.nan
                    salt = np.zeros(N) * np.nan
                    ocean_irr_sol = (np.zeros(N) * np.nan, np.zeros(N) * np.nan, np.zeros(N) * np.nan, np.zeros(N) * np.nan)

                else: 
                    ## Storing the surface chla as the insitu comparison. 
                    ## Calculating the irradiance
                    phy = OI.Phy(z, chla, ESD(phy_type), abscat(lam, phy_type, C2chla='default')[0], abscat(lam, phy_type, C2chla='default')[1])
                    #cdom = OI.CDOM(z, salt, lam)
                    ocean_irr_sol = OI.ocean_irradiance(PI, 
                                                  z[0], 
                                                  abscat(lam, 'water'), 
                                                  phy=phy, 
                                                  N=N)

                                                                
                irr_arr[:,k,0] = ocean_irr_sol[0]
                irr_arr[:,k,1] = ocean_irr_sol[1]
                irr_arr[:,k,2] = ocean_irr_sol[2]
                irr_arr[:,k,3] = ocean_irr_sol[3] 
       
            ## save irr_arr into dict 
            irr_field[lam] = np.copy(irr_arr) 

        ## save to pickle
        pickle.dump(irr_field, open(save_path, "wb"))
        print('Python calculation complete and saved')

    ## If the save file does exist then load it. 
    elif os.path.exists(save_path) == True:
        irr_field = pickle.load(open(save_path,'rb'))
        print('Yay, file loaded :)')

    
    ## Now to caluclate chla. 
    ## Making the Eu at surface dict for ocean color function.
    Rrs = {}
    for lam in wavelengths:
        Eu_surf = np.zeros(N_cst)
        for k,Eu0 in enumerate(irr_field[lam][-1, :, 2]): 
            if Eu0<0 or Eu0>1: 
                Eu0 = np.nan
            Eu_surf[k] = Eu0
        ## Calculating the Rrs.
        Rrs[lam] = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf)

    ## calculating the chla 
    chla_irr = OIR.OCx_alg(Rrs[443], Rrs[490], Rrs[510], Rrs[560]) 
    ## chla data from calcofi 
    chla_dat = np.zeros(N_cst)
    zbot_dat = np.zeros(N_cst)
    lon_dat = cal_cast_dat_bnd['Lon_Dec']
    lat_dat = cal_cast_dat_bnd['Lat_Dec']
    for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
        z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
        if len(z) == 0: 
            z = np.zeros(N) * np.nan
            chla = np.zeros(N) * np.nan
        ## surface is the insitu.
        chla_dat[k] = chla[-1]
        zbot_dat[k] = z[0] 

#!#! Please move Plotting to sperate function to loop over the species.        
 
    ## Plotting the comparison
    if plot:
        fig, ax =plt.subplots()
        ax = PC.Plot_Comparison(ax, chla_dat, chla_irr, 'Chla In Situ to Chla Irr', phy_type, 'In Situ', 'Model', xlim=50, ylim=50)
        ax.legend(title='species')
        fig.show()
    
        ## Plotting the comparison as a scatter on map. 
   #     PF.Plot_Scatter(chla_dat - chla_irr, lat_dat, lon_dat, 'Chla Bias', 'Chla_dat - Chla_irr', vmin=None, vmax=None, fig=None)
    
        ## Plotting the frequency of diff chla vals
        #fig, ax = plt.subplots()
        #N_bins = 1000
        #ax, bin_edges = PC.Plot_Frequency(ax, chla_dat, N_bins, 'Chla_dat') 
        #ax, bin_edges = PC.Plot_Frequency(ax, chla_irr, N_bins, 'Chla_irr', bin_edges=bin_edges) 
        #ax.legend()
        #ax.set_title(f'Frequency Distribtion out of {len(chla_irr)}')
        #ax.set_xlabel('Chla Value [mg chla m^-3]')
        #fig.show()
    
    return chla_dat, irr_field, chla_irr, Rrs


def Loop_Species_Irr_Comp_Cal(PI, save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, species): 
    """
    """

    fig, ax = plt.subplots()

    for k, phy_type in enumerate(species):

        chla_dat, irr_field, chla_irr, Rrs_443, Rrs_551 = Run_Irr_Comp_Insitu(PI, save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, phy_type, plot=False)
        ax = PC.Plot_Comparison(ax, chla_dat, chla_irr, 'Chla In Situ to Chla Irr', phy_type, 'In Situ', 'Model', xlim=50, ylim=50)

    ax.legend()
    fig.show()
  
    return


def Run_Cal_Comp_Viirs(year_min, cal_cast_dat, cal_bot_dat, cci_url, save_dir, save_head, PI, N, wavelengths, phy_type, plot=False):
    """
    The main goal of this function is to plot a one to one comparison of the satellite derived 
    chla values and that from cal cofi. 
    """
    
    save_path = f'{save_dir}/{save_head}'

    ## Bounded data set in the desired timeline.
    cal_cast_dat_bnd = cal_cast_dat[cal_cast_dat['Year'] > year_min]
    ## The number of casts to be calculated. 
    N_cst = len(cal_cast_dat_bnd)

    ## calcofi lon and lat for bounded points.
    lon = cal_cast_dat_bnd['Lon_Dec'].to_numpy()
    lat = cal_cast_dat_bnd['Lat_Dec'].to_numpy()
    year = cal_cast_dat_bnd['Year'].to_numpy()
    julian_day = cal_cast_dat_bnd['Julian_Day'].to_numpy()

    if os.path.exists(f'{save_path}_cal_chla.p') == False:

    
        cal_chla = np.zeros(N_cst)

        ## The Plymouth Marine Lab Data. 
        cci_chla = np.zeros(N_cst)
        cci_Rrs443 = np.zeros(N_cst)
        cci_Rrs490 = np.zeros(N_cst)
        cci_Rrs510 = np.zeros(N_cst)
        cci_Rrs560 = np.zeros(N_cst)
        cci_lat = np.zeros(N_cst)
        cci_lon = np.zeros(N_cst)
        ## Setting the cal cofi limits for the downloaded domain.
        year_lims = [year[0], year[-1]]
        julian_date_lims = [julian_day[0], julian_day[-1]]
        lat_lims = [lat[0], lon[-1]]
        lon_lims = [lon[0], lon[-1]]
        ## Downloading the sub plymouth data set correspoding to the cal cofi domain.
        cci_ds = cci_oc.Get_PML_OC_Data_Set(cci_url, year_lims, julian_date_lims, lat_lims, lon_lims) 
        
        
        ## Now loop over the casts and calculate the irradiance chla each time.  
        for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
            
             z, c_chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt, incld_salt=False)     
             ## checking for zero sized return
             if len(z) == 0: 
                 z = np.zeros(N) * np.nan
                 c_chla = np.zeros(N) * np.nan

             cal_chla[k] = c_chla[-1]
             
             ## Getting the PML data.
             p_Rrs, p_chla, p_lat, p_lon = cci_oc.Get_Point_PML_Dataset(cci_ds, year[k], julian_day[k], lat[k], lon[k]) 
             cci_chla[k] = p_chla
             cci_Rrs443[k] = p_Rrs[443]
             cci_Rrs490[k] = p_Rrs[490]
             cci_Rrs510[k] = p_Rrs[510]
             cci_Rrs560[k] = p_Rrs[560]
             cci_lat[k] = p_lat
             cci_lon[k] = p_lon


        pickle.dump(cal_chla, open(f'{save_path}_cal_chla.p', "wb"))
        pickle.dump(cci_chla, open(f'{save_path}_cci_chla.p', 'wb'))
        pickle.dump(cci_Rrs443, open(f'{save_path}_cci_Rrs443.p', 'wb'))
        pickle.dump(cci_Rrs490, open(f'{save_path}_cci_Rrs490.p', 'wb'))
        pickle.dump(cci_Rrs510, open(f'{save_path}_cci_Rrs510.p', 'wb'))
        pickle.dump(cci_Rrs560, open(f'{save_path}_cci_Rrs560.p', 'wb'))
        pickle.dump(cci_lat, open(f'{save_path}_cci_lat.p', 'wb'))
        pickle.dump(cci_lon, open(f'{save_path}_cci_lon.p', 'wb'))


    elif os.path.exists(f'{save_path}_cal_chla.p') == True:
        
        ## Loading data from pickle.
        cal_chla = pickle.load(open(f'{save_path}_cal_chla.p', "rb"))
        cci_chla = pickle.load(open(f'{save_path}_cci_chla.p', 'rb'))
        cci_Rrs443 = pickle.load(open(f'{save_path}_cci_Rrs443.p', 'rb'))
        cci_Rrs490 = pickle.load(open(f'{save_path}_cci_Rrs490.p', 'rb'))
        cci_Rrs510 = pickle.load(open(f'{save_path}_cci_Rrs510.p', 'rb'))
        cci_Rrs560 = pickle.load(open(f'{save_path}_cci_Rrs560.p', 'rb'))
        cci_lat = pickle.load(open(f'{save_path}_cci_lat.p', 'rb'))
        cci_lon = pickle.load(open(f'{save_path}_cci_lon.p', 'rb'))

    ## The irradiance calculation of Rrs and chla
    cal_chla, irr_field, irr_chla, irr_Rrs = Run_Irr_Comp_Insitu(PI, 
                                                                 save_dir,
                                                                 f'{save_head}_{year_min}',
                                                                 wavelengths,
                                                                 N,
                                                                 year_min,
                                                                 cal_cast_dat,
                                                                 cal_bot_dat,
                                                                 phy_type)
 

    if plot:
        ## Plotting one to one comparison.
        #fig, ax = plt.subplots()
        #PC.Plot_Comparison(ax, cal_chla, viirs_chla, 'Comparison to Cal Chla', 'VIIRS', 'Calcofi Chla', 'VIIRS Chla') 
        #PC.Plot_Comparison(ax, cal_chla, irr_chla, 'Comparison to Cal Chla', 'Irradiance' , 'Calcofi Chla', 'Chla') 
        #ax.legend()
        #fig.show() 
    
        ## Rrs Comparison
        #fig, ax = plt.subplots()
        #PC.Plot_Comparison(ax, viirs_Rrs443, irr_Rrs443, 'Comparison of VIIRS to Irradiance Rrs', '443', 'VIIRS', 'Irr') 
        #PC.Plot_Comparison(ax, viirs_Rrs551, irr_Rrs551, 'Comparison of VIIRS to Irradiance Rrs', '551', 'VIIRS', 'Irr') 
        #ax.legend()  
        #fig.show() 
    
    
        ## Plotting Comparison of chla and location
        fig = plt.figure()
        
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
        #ax.gridlines()
        ## size of scatter points.
        s = 15
        vmin = 0
        vmax = 3
        cbar_shrink = 1
        cbar_label = 'chla [mg chla m^-3]'
    
        ## Plotting the cal cofi
        im = ax.scatter(lon, lat, c=cal_chla, s=s, cmap='nipy_spectral', 
                         transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin, label='calcofi' )
        ## The viirs data.
        im = ax.scatter(cci_lon, cci_lat, c=cci_chla, s=s, cmap='nipy_spectral', 
                         transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin, label='CCI', marker='v')
        ## plotting the lins between corresponding points. 
        for k in range(N_cst):
            plt.plot([lon[k], cci_lon[k]], [lat[k], cci_lat[k]],linewidth=.5, c ='k')
    
        fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = cbar_label)
        ax.set_title('CalCOFI Profile Positions')  
        ylims = ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
        ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels=False
        gl.right_labels=False
        #ax.legend(title=r'mg Chl-a $\mathrm{m}^-1$,)

        fig.show()
 
    
    for k in range(N_cst): 
        dist = geopy.distance.distance((lat[k], lon[k]), (cci_lat[k], cci_lon[k])).kilometers
        if dist > 3:
            cci_chla[k] = np.NaN
            cci_Rrs443[k] = np.NaN
            cci_Rrs490[k] = np.NaN
            cci_Rrs510[k] = np.NaN
            cci_Rrs560[k] = np.NaN
            cci_lat[k] = np.NaN
            cci_lon[k] = np.NaN

    if plot:
        ## Plotting one to one comparison.
        #fig, ax = plt.subplots()
        #PC.Plot_Comparison(ax, cal_chla, viirs_chla, 'Comparison to Cal Chla', 'VIIRS', 'Calcofi Chla', 'VIIRS Chla') 
        #PC.Plot_Comparison(ax, cal_chla, irr_chla, 'Comparison to Cal Chla', 'Irradiance' , 'Calcofi Chla', 'Chla') 
        #ax.legend()
        #fig.show() 
    
        ## Rrs Comparison
        #fig, ax = plt.subplots()
        #PC.Plot_Comparison(ax, viirs_Rrs443, irr_Rrs443, 'Comparison of VIIRS to Irradiance Rrs', '443', 'VIIRS', 'Irr') 
        #PC.Plot_Comparison(ax, viirs_Rrs551, irr_Rrs551, 'Comparison of VIIRS to Irradiance Rrs', '551', 'VIIRS', 'Irr') 
        #ax.legend()  
        #fig.show() 
    
    
        ## Plotting Comparison of chla and location
        fig = plt.figure()
        
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
        #ax.gridlines()
        ## size of scatter points.
        s = 15
        vmin = 0
        vmax = 3
        cbar_shrink = 1
        cbar_label = 'chla [mg chla m^-3]'
    
        ## Plotting the cal cofi
        im = ax.scatter(lon, lat, c=cal_chla, s=s, cmap='nipy_spectral', 
                         transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin, label='calcofi' )
        ## The viirs data.
        im = ax.scatter(cci_lon, cci_lat, c=cci_chla, s=s, cmap='nipy_spectral', 
                         transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin, label='CCI', marker='v')
        ## plotting the lins between corresponding points. 
        for k in range(N_cst):
            plt.plot([lon[k], cci_lon[k]], [lat[k], cci_lat[k]],linewidth=.5, c ='k')
    
        fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = cbar_label)
        ax.set_title('CalCOFI Profile Positions After Location Filter')  
        ylims = ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
        ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels=False
        gl.right_labels=False
        #ax.legend(title=r'mg Chl-a $\mathrm{m}^-1$,)

        fig.show()
 

    cci_Rrs = {}
    cci_Rrs[443] = cci_Rrs443
    cci_Rrs[490] = cci_Rrs490
    cci_Rrs[510] = cci_Rrs510
    cci_Rrs[560] = cci_Rrs560

    return cal_chla,  cci_chla, cci_Rrs, irr_chla, irr_Rrs


def Loop_Species_Viirs_Comp_Cal(year_min, cal_cast_dat, cal_bot_dat, cci_url, save_dir, save_head, PI, N, wavelengths, species): 
    """
    """
    ## The limits for the Rrs plots, in order [xlim, ylim].
    Rrs_ax_lims = [[.014, .014], [.014, .014], [.013, .013], [.0225, .0225], [.014, .014]]

    ncols = len(species)
    fig, axes = plt.subplots(nrows=1, ncols=ncols)
    ## The zoomed in chla plot with different species.
    fig_z, ax_z = plt.subplots()
    axes_list = axes.flatten()
    #chla_ax = axes_list[0]
    #rrs_axes = axes_list[1:]

    for k, phy_type in enumerate(species):
        cal_chla, cci_chla, cci_Rrs, irr_chla, irr_Rrs = Run_Cal_Comp_Viirs(year_min, cal_cast_dat, cal_bot_dat, cci_url, save_dir, save_head, PI, N, wavelengths, phy_type, plot=False)
        
        ## Plotting one to one comparison.
        #chla_ax = PC.Plot_Comparison(chla_ax, cal_chla, irr_chla, 'Comparison to Cal Chla', f'Irr {phy_type}' , 'Calcofi Chla', 'Chla') 
        ax_z = PC.Plot_Comparison(ax_z, cal_chla, irr_chla, 'Irradiance Model Comparison to CalCOFI', f'Irr {phy_type}' , r'CalCOFI Chl-a [mg Chl-a $\mathrm{m}^{-3}$]', '[mg Chl-a $\mathrm{m}^{-3}$]') 
        #ax.legend()
    
        ## Rrs Comparison
        rrs_ax = axes_list[k]

        if k == 0: 
            ylabel = 'Irradiance Model'
        else: 
            ylabel = None
            
        PC.Plot_Comparison(rrs_ax, 
                           cci_Rrs[443], 
                           irr_Rrs[443], 
                           f'{phy_type} Rrs ' + r'[$\mathrm{sr}^{-1}$]', '443 [nm]', 
                           'CCI', 
                           ylabel, 
                           xlim = Rrs_ax_lims[k][0], 
                           ylim=Rrs_ax_lims[k][1], 
                           color= W2RGB.wavelength_to_rgb(443)) 

        PC.Plot_Comparison(rrs_ax, 
                           cci_Rrs[490], 
                           irr_Rrs[490], 
                           f'{phy_type} Rrs ' + r'[$\mathrm{sr}^{-1}$]', 
                           '490 [nm]', 'CCI', 
                           ylabel, 
                           xlim = Rrs_ax_lims[k][0], 
                           ylim=Rrs_ax_lims[k][1], 
                           color=W2RGB.wavelength_to_rgb(490)) 
        PC.Plot_Comparison(rrs_ax, 
                           cci_Rrs[510], 
                           irr_Rrs[510], 
                           f'{phy_type} Rrs ' + r'[$\mathrm{sr}^{-1}$]', 
                           '510 [nm]', 
                           'CCI', 
                           ylabel, 
                           xlim = Rrs_ax_lims[k][0], 
                           ylim=Rrs_ax_lims[k][1], 
                           color =W2RGB.wavelength_to_rgb(510))

        PC.Plot_Comparison(rrs_ax, 
                           cci_Rrs[560], 
                           irr_Rrs[560], 
                           f'{phy_type} Rrs ' + r'[$\mathrm{sr}^{-1}$]', 
                           '551 [nm]', 
                           'CCI', 
                           ylabel, 
                           xlim = Rrs_ax_lims[k][0], 
                           ylim=Rrs_ax_lims[k][1], 
                           color = W2RGB.wavelength_to_rgb(560)) 
        rrs_ax.legend()  
 
    #chla_ax = PC.Plot_Comparison(chla_ax, cal_chla, viirs_chla, 'Comparison to Cal Chla', 'VIIRS (from Rrs/OCx)', 'Calcofi Chla', 'Chla', marker='v', color='grey') 
    #chla_ax = PC.Plot_Comparison(chla_ax, cal_chla, cci_chla, 'Comparison to Cal Chla', 'CCI Chla', 'Calcofi Chla', 'Chla', marker='v', color='black') 
    #ax_z = PC.Plot_Comparison(ax_z, cal_chla, viirs_chla, 'Comparison to Cal Chla', 'VIIRS (from Rrs/OCx)', 'Calcofi Chla', 'Chla', marker='v', color='grey') 
    ax_z = PC.Plot_Comparison(ax_z, cal_chla, cci_chla, 'Irradiance Model Comparison to CalCOFI', f'CCI Chl-a' , r'CalCOFI Chl-a [mg Chl-a $\mathrm{m}^{-3}$]', '[mg Chl-a $\mathrm{m}^{-3}$]', marker='v', color='black') 
    ax_z.legend()
    #chla_ax.legend()
    ## Adding the CCI data info to the figure
    #text(0, .1, 'CCI = ESA CCI Ocean Colour Product Level-3 v5.0 Daily', figure=fig)
    fig_z.show()
    fig.show()

    ## This plot shows a comparison between the CCI calculated chlor and the CCI chlor calculated using
    ## the OCx algorithim from CCI Rrs values. 
    fig, ax = plt.subplots()
    cci_OCx_chla = OIR.OCx_alg(cci_Rrs[443], cci_Rrs[490], cci_Rrs[510], cci_Rrs[560])
    ax = PC.Plot_Comparison(ax, cci_chla, cci_OCx_chla, '', 'CCI', r'CCI Chl-a [mg Chl-a $\mathrm{m}^{-3}$]', r'OC4 Chl-a Using CCI $\mathrm{R_{rs}}$ [mg Chl-a $\mathrm{m}^{-3}$]', marker='v') 
    fig.show()



    return



def Comp_Nomad_Viirs_Irr_Cal(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type):
    """
    This file uses the 13 points from the NOMAD validation data set and compares chla to 
    the satellite value of 
    """

    ## Calculating the irradiance and cal cofi insitu
    irr_rrs443, irr_rrs551, irr_chla, cal_chla = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla='default')


    ## The NOMAD Sattelite and insitu. 
    nomad_insitu_chla = chla_val_cal_dat['insitu_chlor_a']
    nomad_sat_chla = chla_val_cal_dat['aqua_chlor_a']
    nomad_sat_rrs_443 = chla_val_cal_dat['aqua_rrs443'] 
    nomad_sat_rrs_555 = chla_val_cal_dat['aqua_rrs555'] 
    nomad_lon = chla_val_cal_dat['longitude'].to_numpy()
    nomad_lat = chla_val_cal_dat['latitude'].to_numpy()
    nomad_date_time = chla_val_cal_dat['date_time']

    ## The chla and rrs from viirs
    ## VIIRS went only started in 2012.
    ## The chla_val_cal_dat in around 2003
#    viirs_chla = np.zeros(len(chla_val_cal_dat))
#    viirs_rrs443 = np.zeros(len(chla_val_cal_dat))
#    viirs_rrs551 = np.zeros(len(chla_val_cal_dat))
#    
#    for k, date_time in enumerate(chla_val_cal_dat['date_time']):
#        date, time = date_time.split()
#        year = date.split('-')[0]
#        julian_day =  Date_to_Julian_Date(date)   
#
#        rrs443, rrs551, chla = viirs_calcofi_val.Viirs_Rrs_Chla(julian_day, year, nomad_lat[k],nomad_lon[k])
#        
#        viirs_chla[k] = chla
#        viirs_rrs443[k] = rrs443
#        viirs_rrs551[k] = rrs551

    ## Plotting
    ## --------

    ## Plotting the chla's against nomad in situ chla.
    fig, ax = plt.subplots()

    title = 'Chla Comparisons'
    xlabel = 'Chla NOMAD In Situ [mg chla m^-3]'
    ylabel = 'Chla [mg chla m^-3]'
    ax = PC.Plot_Comparison(ax, nomad_insitu_chla, nomad_sat_chla, title, 'NOMAD Satellite', xlabel, ylabel)
    ax = PC.Plot_Comparison(ax, nomad_insitu_chla, cal_chla , title, 'Calcofi In Situ', xlabel, ylabel)
    ax = PC.Plot_Comparison(ax, nomad_insitu_chla, irr_chla , title, 'Irr Model', xlabel, ylabel)
#    ax1 = PC.Plot_Comparison(ax1, nomad_insitu_chla, viirs_chla, title, 'VIIRS', xlabel, ylabel)
    ax.legend()

    ## Plotting the rrs's against the nomad rrs.
    ## 443
#    title = 'Rrs 443 Comparison'
#    xlabel = 'Rrs 443 NOMAD Satellite [sr^-1]'
#    ylabel = 'Rrs 443 [sr^-1]'
#    ax2 = PC.Plot_Comparison(ax2, nomad_sat_rrs_443, irr_rrs443, title, 'Irr Model', xlabel, ylabel)
#    ax2 = PC.Plot_Comparison(ax2, nomad_sat_rrs_443, viirs_rrs443, title, 'VIIRS', xlabel, ylabel)
#    ax2.legend()

    ## 551
#    title = 'Rrs 551 Comparison'
#    xlabel = 'Rrs 551 NOMAD Satellite [sr^-1]'
#    ylabel = 'Rrs 551 [sr^-1]'
#    ax3 = PC.Plot_Comparison(ax3, nomad_sat_rrs_555, irr_rrs551, title, 'Irr Model', xlabel, ylabel)
#    ax3 = PC.Plot_Comparison(ax3, nomad_sat_rrs_555, viirs_rrs551, title, 'VIIRS', xlabel, ylabel)
#    ax3.legend()
 
    fig.show()

   
    return 
    

def Count_Cal_Chla_Samp_Pts(year_min, cal_cast_dat, cal_bot_dat):
    """
    This file counts the number of sampling point used in each calcofi cast and plots 
    it against the value of chla at the surface.
    """

    ## Bounded data set in the desired timeline.
    cal_cast_dat_bnd = cal_cast_dat[cal_cast_dat['Year'] > year_min]
    ## The number of casts to be calculated. 
    N_cst = len(cal_cast_dat_bnd)

    ## The number of points sampled fopr each cast
    cst_samps = np.zeros(N_cst)
    cst_zbot = np.zeros(N_cst)
    cst_chla = np.zeros(N_cst)

    for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
            
         z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt, incld_salt=False)     
         ## checking for zero sized return
         if len(z) == 0: 
             cst_samps[k] = 0
             cst_zbot[k] = np.nan 
             cst_chla[k] = np.nan
         else: 
             cst_samps[k] = len(chla)
             cst_zbot[k] = z[0]
             cst_chla[k] = chla[-1]
       
    ## Plotting the results
    fig, axes =  plt.subplots(nrows = 1, ncols = 2)

    ## Plotting the number of sampling points against chla.
    axes[0].plot(cst_chla, cst_samps, 'o',fillstyle='none')
    axes[0].set_xlabel('Chla')
    axes[0].set_ylabel('Number of Sampling Points in Cast')
    axes[0].set_title('Number of Sampling Points')

    ## Plotting the bottom depth of cast. 
    axes[1].plot(cst_chla, cst_zbot, 'o',fillstyle='none')
    axes[1].set_xlabel('Chla')
    axes[1].set_ylabel('Zbot [m]')
    axes[1].set_title('Bottom depth of the Cast')
    
    fig.show()

    return 


if __name__ == '__main__': 

    import argparse 
    parser = argparse.ArgumentParser(description='Ocean Irradiance comparisons to CalCofi cruise data')
    parser.add_argument('save_dir', help='The directory of the save file.')
    parser.add_argument('save_file_head', help='The name of the pickle file to be saved. Shoud not include the .p')
    args = parser.parse_args()

    ## The species that have coefficients.
    #species = ['Diat', 'Cocco']
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']
    #species = ['HLPro', 'Cocco', 'Diat']
 
    ## The spread of different C2chla ratios. 
    C2chla_vals = np.arange(50,200, 25)

    ## The data from files. 
    chla_val_dat, cal_cast_dat, cal_bot_dat = Get_Data('chlor_a_validation.csv', 'cal_casts.csv', 'cal_bottles.csv')
    
    ## The url of the plymouth data
    ## From the website : https://www.oceancolour.org/thredds/catalog-cci.html?dataset=CCI_ALL-v5.0-DAILY
    ## time is in units of days since 1970-01-01 00:00:00
    #ply_dat_url = 'https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v5.0-DAILY?lat[0:1:0],lon[0:1:0],time[0:1:0],Rrs_443[0:1:0][0:1:0][0:1:0],Rrs_560[0:1:0][0:1:0][0:1:0],chlor_a[0:1:0][0:1:0][0:1:0]' 

    cci_url = 'https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v5.0-DAILY?lat[0:1:0],lon[0:1:0],time[0:1:0],Rrs_443[0:1:0][0:1:0][0:1:0],Rrs_490[0:1:0][0:1:0][0:1:0],Rrs_510[0:1:0][0:1:0][0:1:0],Rrs_560[0:1:0][0:1:0][0:1:0],chlor_a[0:1:0][0:1:0][0:1:0]' 
#    'https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v5.0-DAILY?lat[0:1:4319],lon[0:1:8639],time[0:1:8858],Rrs_443[0:1:0][0:1:0][0:1:0],Rrs_490[0:1:0][0:1:0][0:1:0],Rrs_510[0:1:0][0:1:0][0:1:0],Rrs_560[0:1:0][0:1:0][0:1:0]'
    ## getting the data with cal cruises only
    chla_val_cal_dat = chla_val_dat[chla_val_dat['cruise'].str.contains('cal', regex = False) ]
 
    ## Angle of the azimuth. Only for baird. 
    theta_air=.55

    ## Running 
    #Run_and_Plot_Comparison(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, species, C2chla_vals)

    ## The number of vertical layers in irr grid. 
    N = 500

    ## The oldest.
    ## should be around 5000 casts. 
    year_min = 2012
 
    ## Wavelengths
    wavelengths = [443, 490, 510, 560]

    ## The species of phytoplankton
    phy_type = 'Diat'

    ## The save file designated by the desired min year. 
    save_file = f'{args.save_file_head}_{year_min}'
  
    ## Param Init object 
    PI = Param_Init()

    ## Running comparison of insitu to irr surface chla for many cal casts. 
    #Run_Irr_Comp_Insitu(PI, args.save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, 'Diat', plot=True)
    #Loop_Species_Irr_Comp_Cal(PI, args.save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, species)

    save_path = f'{args.save_dir}/{args.save_file_head}'
    phy_type = 'Diat'
    ## Runnning the comparison of calcofi to viirs
    Run_Cal_Comp_Viirs(year_min, cal_cast_dat, cal_bot_dat, cci_url, args.save_dir, args.save_file_head, PI, N, wavelengths, phy_type, plot=True) 
#    Loop_Species_Viirs_Comp_Cal(year_min, cal_cast_dat, cal_bot_dat, cci_url, args.save_dir, args.save_file_head, PI, N, wavelengths, species)


    ## Running the comparison of viirs, calcofi, irr, and nomad
#    Comp_Nomad_Viirs_Irr_Cal(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type)
   
    ## Plotting the sample points of th casts
    #Count_Cal_Chla_Samp_Pts(year_min, cal_cast_dat, cal_bot_dat)
