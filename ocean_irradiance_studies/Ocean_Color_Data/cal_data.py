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
## User made mod
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
import ocean_irradiance_visualization.Plot_Field as PF


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


def Get_Cast_Depth_Chla(cal_bot_dat, cast_count): 
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
    z = z[~np.isnan(chla)]
    chla = chla[~np.isnan(chla)]

    return z, chla

 
def Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=None): 
    """
    This function loops over the cruise data and calculates the irradiances and chla at the surface for each. 
    """ 



    ## using my params
    PI = Param_Init()

    Eu_surf = np.zeros(len(chla_val_cal_dat['cruise']))
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
            z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
            ## Calculating the irradiance
            phy = OI.Phy(z, chla, abscat(lam, phy_type, C2chla=C2chla)[0], abscat(lam, phy_type, C2chla=C2chla)[1])
            irr_out = OI.ocean_irradiance_shoot_up(
                                                   z[0],
                                                   PI.Ed0,
                                                   PI.Es0,
                                                   PI.Euh,
                                                   abscat(lam, 'water'),
                                                   PI.coefficients,
                                                   phy=phy,
                                                   N=100,
                                                   pt1_perc_zbot=True
                                                   )
            Eu_surf[k] =  irr_out[2][-1]
        Eu_surf_dict[lam] = np.copy(Eu_surf)
            
         
    ## Calculating the rrs for each wavelength and the chla 
    rrs_443 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf_dict[443])
    rrs_551 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf_dict[551])

    ## Chla 
    irr_chla = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)

    return rrs_443, rrs_551, irr_chla        
            

def Plot_Comparison(ax, x, y, title, label, xlabel, ylabel, xlim=None, ylim=None): 
    """
    Plots the given values on a given axes
    """

    ax.plot(x, y,'o', fillstyle='none', label=label, markersize=5)
    ax.plot(x, x, 'k')
    if xlim == None: 
        xlim = x.max()
    if ylim == None: 
        ylim = y.max()
#    ax.set_xlim([0, xlim])
#    ax.set_ylim([0, ylim])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    return ax 


def Run_and_Plot_Comparison(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, species, C2chla_vals):
    """
    Plots a scatter plot comparing the values
    """ 

    ## The reference data. 
    chla_dat = chla_val_cal_dat['insitu_chlor_a']
    chla_sat = chla_val_cal_dat['aqua_chlor_a']
    rrs_443_dat = chla_val_cal_dat['aqua_rrs443'] 
    rrs_555_dat = chla_val_cal_dat['aqua_rrs555'] 
    rrs_ratio_dat =  rrs_443_dat / rrs_555_dat 

    fig, axes = plt.subplots(3, 2)
    
    ## Plot C2chla as LC2chla ratio with diff species. 
    LC2chla = 100
    for k, phy_type in enumerate(species):
        rrs_443, rrs_551, irr_chla = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=LC2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## rrs slope 
        rrs_m_443 = rrs_443 / rrs_443_dat
        rrs_m_551 = rrs_551 / rrs_555_dat
        ## chla comparison
        Plot_Comparison(axes[0,0], chla_dat, irr_chla, 'LC2chla Ratio Varied Species Chla', phy_type, 'Insitu', 'Model') 
        Plot_Comparison(axes[0,1], rrs_ratio_dat, rrs_ratio, 'LC2chla Ratio Varied Species RRS Ratio', phy_type, 'Insitu', 'Model') 
    ## plotting satelite 
    Plot_Comparison(axes[0,0], chla_dat, chla_sat, 'LC2chla Ratio Varied Species Chla', 'Sattelite', 'Insitu', 'Model') 
    axes[0,0].legend( title='species')
    axes[0,0].grid()
    axes[0,1].grid()

    ## Plot C2chla as SC2chla ratio with diff species. 
    SC2chla = 50
    for k, phy_type in enumerate(species):
        rrs_443, rrs_551, irr_chla = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=SC2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        Plot_Comparison(axes[1,0], chla_dat, irr_chla, 'SC2chla Ratio Varied Species Chla', phy_type, 'Insitu', 'Model') 
        Plot_Comparison(axes[1,1], rrs_ratio_dat, rrs_ratio, 'SC2chla Ratio Varied Species RRS Ratio', phy_type, 'Insitu', 'Model') 
    Plot_Comparison(axes[1,0], chla_dat, chla_sat, 'LC2chla Ratio Varied Species Chla', 'Sattelite', 'Insitu', 'Model') 
    axes[1,0].legend(title='species')
    axes[1,0].grid()
    axes[1,1].grid()

    ## Plot generic species with different C2chla. 
    phy_type = 'Generic'
    for k, C2chla in enumerate(C2chla_vals):
        rrs_443, rrs_551, irr_chla = Loop_Cal_Cruise(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, C2chla=C2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        Plot_Comparison(axes[2,0], chla_dat, irr_chla, 'Generic Species Varied C2chla Chla', C2chla, 'Insitu', 'Model') 
        Plot_Comparison(axes[2,1], rrs_ratio_dat, rrs_ratio, 'Generic Species Varied C2chla RRS Ratio', C2chla, 'Insitu', 'Model') 
    Plot_Comparison(axes[2,0], chla_dat, chla_sat, 'LC2chla Ratio Varied Species Chla', 'Sattelite', 'Insitu', 'Model') 
    axes[2,0].legend(title='C2chla')
    axes[2,0].grid()
    axes[2,1].grid()

    #plt.tight_layout()
    fig.show()


    return rrs_m_443, rrs_m_551 
 

def Run_Irr(PI, save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, phy_type):
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
    
    ## The location of the save file.
    save_path = f'{save_dir}/{save_file}'
 
    ## Bounded data set in the desired timeline.
    cal_cast_dat_bnd = cal_cast_dat[cal_cast_dat['Year'] > year_min]
    ## The number of casts to be calculated. 
    N_cst = len(cal_cast_dat_bnd)

    ## The ditionary of irr fields. 
    irr_field = {}

    ## Checking if the save file exists.
    ## If it doesn't then do the calculation.
    if os.path.exists(save_path) == False:

        for lam in wavelengths: 
            ## The array in which to store all the irradiance solutions into. 
            irr_arr = np.zeros((N,N_cst, 4))
         
            ## Now loop over the casts and calculate the irradiance chla each time.  
            for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
                print(f'{k}/{N_cst}')
                ## Getting the chla and depth profile.
                z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
                ## checking for zero sized return
                if len(z) == 0: 
                    z = np.zeros(N) * np.nan
                    chla = np.zeros(N) * np.nan
                ## Storing the surface chla as the insitu comparison.

                ## Calculating the irradiance
                phy = OI.Phy(z, chla, abscat(lam, phy_type)[0], abscat(lam, phy_type)[1])
                ocean_irr_sol = OI.ocean_irradiance_shoot_up(
                                                             z[0],
                                                             PI.Ed0,
                                                             PI.Es0,
                                                             PI.Euh,
                                                             abscat(lam, 'water'),
                                                             PI.coefficients,
                                                             phy=phy,
                                                             N=N,
                                                             pt1_perc_zbot=True
                                                             )
                ## Ed, Es, Eu, z 
                irr_arr[:,k,0] = ocean_irr_sol[0]
                irr_arr[:,k,1] = ocean_irr_sol[1]
                irr_arr[:,k,2] = ocean_irr_sol[2]
                irr_arr[:,k,3] = ocean_irr_sol[3]
       
            ## save irr_arr into dict 
            irr_field[lam] = irr_arr

        ## save to pickle
        pickle.dump(irr_field, open(save_path, "wb"))
        print('Python calculation complete and saved')

    ## If the save file does exist then load it. 
    elif os.path.exists(save_path) == True:
        irr_field = pickle.load(open(save_path,'rb'))
        print('Yay, file loaded :)')


    ## Now to caluclate chla. 
    ## Making the Eu at surface dict for ocean color function.
    Eu_surf_dict = {}
    Eu_surf = np.zeros(N_cst)
    for lam in wavelengths:
        for k,Eu in enumerate(irr_field[lam][-1, :, 2]): 
            if Eu<0 or Eu>1: 
                Eu_surf[k] = np.nan
            else:
                Eu_surf[k] = Eu
            
        Eu_surf_dict[lam] = np.copy(Eu_surf)
    ## calculating the chla 
    chla_irr = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)

    ## chla data from calcofi 
    chla_dat = np.zeros(N_cst)
    lon_dat = cal_cast_dat_bnd['Lon_Dec']
    lat_dat = cal_cast_dat_bnd['Lat_Dec']
    for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
        z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
        if len(z) == 0: 
            z = np.zeros(N) * np.nan
            chla = np.zeros(N) * np.nan
        ## surface is the insitu.
        chla_dat[k] = np.nanmean(chla[-5:])

    ## Plotting the comparison
    fig, ax = plt.subplots()
    Plot_Comparison(ax, chla_dat, chla_irr, 'Chla Insitu to Chla Irr', phy_type, 'Insitu', 'Model', xlim=None, ylim=None)
    fig.show()

    ## Plotting the comparison as a scatter on map. 
    PF.Plot_Scatter(chla_dat - chla_irr, lat_dat, lon_dat, 'Chla Bias', 'Chla_dat - Chla_irr', vmin=None, vmax=None, fig=None)
    
    return irr_field 
 
        


if __name__ == '__main__': 

    import argparse 
    parser = argparse.ArgumentParser(description='Ocean Irradiance comparisons to CalCofi cruise data')
    parser.add_argument('save_dir', help='The directory of the save file.')
    parser.add_argument('save_file_head', help='The name of the pickle file to be saved. Shoud not include the .p')
    args = parser.parse_args()

    ## The species that have coefficients.
    species = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']
 
    ## The spread of different C2chla ratios. 
    C2chla_vals = np.arange(50,200, 25)

    ## The data from files. 
    chla_val_dat, cal_cast_dat, cal_bot_dat = Get_Data('chlor_a_validation.csv', 'cal_casts.csv', 'cal_bottles.csv')

    ## getting the data with cal cruises only
    chla_val_cal_dat = chla_val_dat[chla_val_dat['cruise'].str.contains('cal', regex = False) ]

    ## Running 
    rrs_m_443, rrs_m_551 = Run_and_Plot_Comparison(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, species, C2chla_vals)

    ## The number of vertical layers in irr grid. 
    N = 100

    ## The oldest that data can be used.
    ## should be around 5000 casts. 
    year_min = 2000
 
    ## Wavelengths
    wavelengths = [443, 551]

    ## The species of phytoplankton
    phy_type = 'Diat'

    ## The save file designated by the desired min year. 
    save_file = f'{args.save_file_head}_{year_min}_{phy_type}.p'
  
    ## Param Init object 
    PI = Param_Init()

    ## Running comparison of insitu to irr surface chla for many cal casts. 
    irr_field = Run_Irr(PI, args.save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, phy_type)  


   
