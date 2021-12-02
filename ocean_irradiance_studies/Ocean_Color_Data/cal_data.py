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
import ocean_irradiance_baird.ocean_irradiance_baird as OIB
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
import ocean_irradiance_visualization.Plot_Field as PF
import ocean_irradiance_visualization.Plot_Comparison as PC


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
            z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
            chla_dat[k] = chla[-1]
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
                                                   N=1000,
                                                   pt1_perc_zbot=True, 
                                                   pt1_perc_phy=False    
                                                   )
            Eu_surf[k] =  irr_out[2][-1]
        Eu_surf_dict[lam] = np.copy(Eu_surf)
            
         
    ## Calculating the rrs for each wavelength and the chla 
    rrs_443 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf_dict[443])
    rrs_551 = OIR.R_RS(PI.Ed0, PI.Es0, Eu_surf_dict[551])

    ## Chla 
    irr_chla = OIR.ocean_color_sol(Eu_surf_dict, PI.Ed0, PI.Es0)

    return rrs_443, rrs_551, irr_chla, chla_dat        


def Loop_Cal_Cruise_Baird(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, theta_air, C2chla=None): 
    """
    This function loops over the cruise data and calculates the Rrs and chla using the 
    Baird model. 
    """ 

    Rrss = np.zeros(len(chla_val_cal_dat['cruise']))
    Rrs_dict = {}
    ## Looping over ocean color wavelengths.
    for lam in [443, 551]:
        ## Then loop over cal cruises. 
        for k, id in enumerate(chla_val_cal_dat['/fields=id']): 
            ## Getting the cast count. 
            cst_cnt = Get_Cst_Cnt(chla_val_cal_dat, cal_cast_dat, id)
            if cst_cnt == None: 
                Rrss[k] = np.nan
                continue 
            ## Getting the chla and depth profile.
            z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
            ## Calculating the irradiance
            phy = OI.Phy(z, chla, abscat(lam, phy_type, C2chla=C2chla)[0], abscat(lam, phy_type, C2chla=C2chla)[1])

            Rrs = OIB.ocean_irradiance_baird(
                                             z[0],
                                             abscat(lam, 'water'),
                                             theta_air,
                                             phy=phy,
                                             N=1000,
                                             )
            Rrss[k] =  Rrs
        Rrs_dict[lam] = np.copy(Rrss)
            
         
    ## Chla 
    chla = OIR.OCx_alg(Rrs_dict[443], Rrs_dict[551])

    return Rrs_dict[443], Rrs_dict[551], chla 
            

def Plot_Comparison(ax, x, y, title, label, xlabel, ylabel, xlim=None, ylim=None): 
    """
    Plots the given values on a given axes
    """

    ax.plot(x, y,'o', fillstyle='none', label=label, markersize=5)
    ax.plot(x, x, 'k')
    if xlim == None: 
        ax.relim()
    elif ylim == None: 
        ax.relim()
    else:
        ax.set_xlim([0, xlim])
        ax.set_ylim([0, ylim])

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    return ax 


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
        if method == 'baird':
            rrs_443, rrs_551,  irr_chla, chla_dat = Loop_Cal_Cruise_Baird(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, theta_air, C2chla=LC2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        print(rrs_443)
        ## chla comparison
        Plot_Comparison(axes[0,0], chla_insitu, irr_chla, 'Chla LC2chla Varied Species', phy_type, None, 'Model') 
        Plot_Comparison(axes[0,1], rrs_ratio_dat, rrs_ratio, 'Rrs Ratio  LC2chla Varied Species', phy_type, None, None) 
        Plot_Comparison(axes2[0,0], rrs_443_dat, rrs_443, 'Rrs 443 LC2chla Varied Species', phy_type, None, 'Model') 
        Plot_Comparison(axes2[0,1], rrs_555_dat, rrs_551, 'Rrs 551 LC2chla Varied Species', phy_type, None, None) 
    ## plotting satelite 
    Plot_Comparison(axes[0,0], chla_insitu, chla_sat, 'Chla LC2chla Varied Species', 'Sattelite', None, 'Model') 
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
        if method == 'baird':
            rrs_443, rrs_551,  irr_chla, chla_dat = Loop_Cal_Cruise_Baird(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, theta_air, C2chla=SC2chla)     
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        Plot_Comparison(axes[1,0], chla_insitu, irr_chla, 'Chla SC2chla Varied Species', phy_type, None, 'Model') 
        Plot_Comparison(axes[1,1], rrs_ratio_dat, rrs_ratio, 'Rrs Ratio SC2chla Varied Species', phy_type, None, None) 
        Plot_Comparison(axes2[1,0], rrs_443_dat, rrs_443, 'Rrs 443 LC2chla Varied Species', phy_type, None, 'Model') 
        Plot_Comparison(axes2[1,1], rrs_555_dat, rrs_551, 'Rrs 551 LC2chla Varied Species', phy_type, None, None) 
    Plot_Comparison(axes[1,0], chla_insitu, chla_sat, 'Chla SC2chla Varied Species', 'Sattelite', None, 'Model') 
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
        if method == 'baird':
            rrs_443, rrs_551,  irr_chla, chla_dat = Loop_Cal_Cruise_Baird(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, phy_type, theta_air, C2chla=C2chla)     
        
        ## Ratio of rrs 
        rrs_ratio = rrs_443 / rrs_551
        ## chla comparison
        Plot_Comparison(axes[2,0], chla_insitu, irr_chla, 'Chla Generic Species Varied C2chla', C2chla, 'Insitu', 'Model') 
        Plot_Comparison(axes[2,1], rrs_ratio_dat, rrs_ratio, 'Rrs Ratio Generic Species Varied C2chla', C2chla, 'Insitu', None) 
        Plot_Comparison(axes2[2,0], rrs_443_dat, rrs_443, 'Rrs 443 LC2chla Varied Species', C2chla, 'Insitu', 'Model') 
        Plot_Comparison(axes2[2,1], rrs_555_dat, rrs_551, 'Rrs 551 LC2chla Varied Species', C2chla, 'Insitu', None) 
    Plot_Comparison(axes[2,0], chla_insitu, chla_sat, 'Chla Generic Species Varied C2chla', 'Sattelite', 'Insitu', 'Model') 
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
    Plot_Comparison(ax, chla_insitu, chla_dat, 'Chla NOMAD Insitu Compared to Cast Insitu', None, 'NOMAD', 'Cast')
    ax.grid()
    fig.show()


    return 
 

def Run_Irr_Comp_Insitu(PI, save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, species):
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
    fig, ax = plt.subplots()
    for phy_type in species:
        ## The location of the save file.
        save_path = f'{save_dir}/{save_file}_{phy_type}.p'

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
                    #if lam == 443: 
                    #    PI.Ed0= .7
                    #    PI.Es0=.3
                    #if lam == 551: 
                    #    PI.Ed0 =.3
                    #    PI.Es0 =.7
                    ocean_irr_sol = OI.ocean_irradiance_shoot_up(
                                                                 z[0],
                                                                 PI.Ed0,
                                                                 PI.Es0,
                                                                 PI.Euh,
                                                                 abscat(lam, 'water'),
                                                                 PI.coefficients,
                                                                 phy=phy,
                                                                 N=N,
                                                                 pt1_perc_zbot=True,
                                                                 pt1_perc_phy=False
                                                                 )
                    ## Ed, Es, Eu, z 
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
        Eu_surf_dict = {}
        f_field = {}
        for lam in wavelengths:
            Eu_surf = np.zeros(N_cst)
            ## calculation of the f paramter from Ocean Optics Webbook
            f_arr = np.zeros((N_cst)) 
            for k,Eu in enumerate(irr_field[lam][-1, :, 2]): 
                if Eu<0 or Eu>1: 
                    Eu_surf[k] = np.nan
                else:
                    Eu_surf[k] = Eu
                  
                ## calculation of f. 
                ## follows the ocean optics webbook normalizing radiances chapter.
                a_wat, b_wat = abscat(lam, 'water')
                b_b_wat = .551 * b_wat
                f = (Eu / (PI.Ed0 + PI.Es0)) * (a_wat/b_b_wat) 
                f_arr[k] = f 
    
            Eu_surf_dict[lam] = np.copy(Eu_surf)
            f_field[lam] = np.copy(f_arr)
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
            chla_dat[k] = chla[-1]
    
        ## Plotting the comparison
        ax = Plot_Comparison(ax, chla_dat, chla_irr, 'Chla Insitu to Chla Irr', phy_type, 'Insitu', 'Model', xlim=50, ylim=50)
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
    
    return irr_field, f_field, chla_dat 
    
    
def Run_Baird_Comp_Insitu(save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, phy_type, theta_air):
    """
    This calculates the Baird chla value for many different casts within a given time line. 
   
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
    """
    
    ## The location of the save file.
    save_path = f'{save_dir}/{save_file}'
 
    ## Bounded data set in the desired timeline.
    cal_cast_dat_bnd = cal_cast_dat[cal_cast_dat['Year'] > year_min]
    ## The number of casts to be calculated. 
    N_cst = len(cal_cast_dat_bnd)

    Rrs_dict = {}
    for lam in wavelengths: 
        ## The array in which to store all the irradiance solutions into. 
        Rrs_arr = np.zeros(N_cst)
         
        ## Now loop over the casts and calculate the irradiance chla each time.  
        for k, cst_cnt in enumerate(cal_cast_dat_bnd['Cst_Cnt'].to_numpy()):
            print(f'{k}/{N_cst}')
            ## Getting the chla and depth profile.
            z, chla = Get_Cast_Depth_Chla(cal_bot_dat, cst_cnt)     
            ## checking for zero sized return
            if len(z) == 0: 
                z = np.zeros(N) * np.nan
                chla = np.zeros(N) * np.nan
            ## Storing the surface chla as the insitu comparis
            ## Calculating the irradiance
            phy = OI.Phy(z, chla, abscat(lam, phy_type)[0], abscat(lam, phy_type)[1])
            Rrs = OIB.ocean_irradiance_baird(
                                             z[0],
                                             abscat(lam, 'water'),
                                             theta_air,
                                             phy=phy,
                                             N=N,
                                             )
            Rrs_arr[k] = Rrs
        ## Save into a dict. 
        Rrs_dict[lam] = np.copy(Rrs_arr)

    ## Now to caluclate chla. 
    chla_irr = OIR.OCx_alg(Rrs_dict[443], Rrs_dict[551])

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
        chla_dat[k] = chla[-1]

    ## Plotting the comparison
    fig, ax = plt.subplots()
    Plot_Comparison(ax, chla_dat, chla_irr, 'Chla Insitu to Chla Irr', phy_type, 'Insitu', 'Model', xlim=None, ylim=None)
    fig.show()

    ## Plotting the comparison as a scatter on map. 
    PF.Plot_Scatter(chla_dat - chla_irr, lat_dat, lon_dat, 'Chla Bias', 'Chla_dat - Chla_irr', vmin=None, vmax=None, fig=None)

    ## Plotting the frequency of diff chla vals
    fig, ax = plt.subplots()
    N_bins = 1000
    ax, bin_edges = PC.Plot_Frequency(ax, chla_dat, N_bins, 'Chla_dat') 
    ax, bin_edges = PC.Plot_Frequency(ax, chla_irr, N_bins, 'Chla_irr', bin_edges=bin_edges) 
    ax.legend()
    ax.set_title(f'Frequency Distribtion out of {len(chla_irr)}')
    ax.set_xlabel('Chla Value [mg chla m^-3]')
    fig.show()

    return  
 
        


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
 
    ## Angle of the azimuth. 
    theta_air=.55

    ## Running 
    #Run_and_Plot_Comparison(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, species, C2chla_vals, method='baird', theta_air=theta_air)
    Run_and_Plot_Comparison(chla_val_cal_dat, cal_cast_dat, cal_bot_dat, species, C2chla_vals)

    ## The number of vertical layers in irr grid. 
    N = 200

    ## The oldest that data can be used.
    ## should be around 5000 casts. 
    year_min = 2000
 
    ## Wavelengths
    wavelengths = [443, 551]

    ## The species of phytoplankton
    phy_type = 'Diat'

    ## The save file designated by the desired min year. 
    save_file = f'{args.save_file_head}_{year_min}'
  
    ## Param Init object 
    PI = Param_Init()

    ## Running comparison of insitu to irr surface chla for many cal casts. 
    #Run_Baird_Comp_Insitu(args.save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, phy_type, theta_air)  
    Run_Irr_Comp_Insitu(PI, args.save_dir, save_file, wavelengths, N, year_min, cal_cast_dat, cal_bot_dat, species)


   
