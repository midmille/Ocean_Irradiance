# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 10:40:56 2021

@author: miles

ABOUT
------

The purpose of this file is to implement the irradiance algorithim with chlor_a 
data from Ocean Color website, https://oceancolor.gsfc.nasa.gov/atbd/chlor_a/ 
using seawifs level 2 satellite-to-in-situ match-up validation data. 
    
"""
## External Modules
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

## User Modules
import Read_CSV_Dat as RCD
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat




def Scatter_Plot_Data_Global(longitude, latitude, val, units, title):
    
    """
    This plots the given value as a scatter plot of color mapped points 
    on a Cartopy module projection. 

    Parameters
    ----------
    longitude : 1-D Array [N]
        Longitude values. 
    latitude : 1-D Array [N]
        latitude values.
    val : 1-D Array [N]
        The value to be plotted. 
    units : String
        The units of the value being plotted.
    title : String
        The title of the plot, the vlaue being plotted. 

    Returns
    -------
    fig : matplotlib Figure instance
        Figure corresponding to this plot.
    ax1 : matplotlib axes instance
        Axes corresponding to this plot.

    """
    
    
    fig = plt.figure(figsize=(8,6))
    cbar_shrink=.5
    
    ## The seawifs chlor
    ax1 = fig.add_subplot(111, projection=ccrs.Robinson())
    ax1.add_feature(cfeature.COASTLINE)
    ax1.gridlines()
    im1 = ax1.scatter(longitude, latitude, c=val, s=5, cmap='nipy_spectral', 
                     transform=ccrs.PlateCarree() )
    fig.colorbar(im1, ax=ax1, shrink=cbar_shrink, label = units)
    ax1.set_title(title)
    
    return fig,ax1


def Mask(long_bound, lat_bound, long, lat, plot=False, ax=None):
    """
    This creates a [0,1] mask for the zone of interest. 

    Parameters
    ----------
    long_bound : 1-D Array [2]
        The longitude boundary for this mask. In increasing order. Array of length 2.
    lat_bound : 1-D Array [2]
        The latitude boundary for this mask. In increasing order. Array of length 2.
    long : 1-D Array [N]
        Longitude array.
    lat : 1-D Array [N]
        latitude array.
    plot : Bool, optional
       Plot mask on a given axes instance. The default is False.
    ax : Bool, optional
        The axes instance to plot mask boundary on. The default is None.

    Returns
    -------
    mask : 1-D Array
        The mask. 1 is in the mask. 0 is out of mask zone. 

    """
    ## Creating a mask for the chl in west coast area.
    
    assert len(long) == len(lat)
    
    ## linewidth 
    lw =2
    
    ## longitude lines
    if plot: 
        ax1 = ax
        
        ax1.plot([long_bound[0],long_bound[0]], lat_bound, 'r', transform=ccrs.PlateCarree(), 
                 linewidth = lw )
        ax1.plot([long_bound[1],long_bound[1]], lat_bound, 'r', transform=ccrs.PlateCarree(), 
                 linewidth = lw )
        ## latitude
        ax1.plot(long_bound, [lat_bound[0],lat_bound[0]], 'r', transform=ccrs.PlateCarree(), 
                 linewidth = lw )
        ax1.plot(long_bound, [lat_bound[1],lat_bound[1]], 'r', transform=ccrs.PlateCarree(), 
                 linewidth = lw, label='MASK' )
        
        ax1.legend()
        
    ## Multiplication of booleans is T*T = 1, T*F = 0, F*F =0
    mask = (lat_bound[0] < lat) * (lat < lat_bound[1]) * (long_bound[0] < long) * (long < long_bound[1])

    return mask


def R_RS_Irr(lam, chlor_dat, use_art_phy=False):
    """
    This calculates the rrs value from an ocean irradiance calculation using a 
    constant phy profile from data. 

    Parameters
    ----------
    lam : Integer
        The given wavelength
    chlor_dat : 1-D Array [N]
        The chlor_a data to be made into a constant phy profile for OI. 

    Returns
    -------
    rrs_irr : 1-D Array [N]
        The irradiance r_rs calculated value to be returned. 

    """
    
    Ed0 = 0.7
    Es0 = 1- Ed0 
    Euh = 0.0
    
    ## The number of levels in profile [same as west coast domain]
    N = 42
    ## The vertical grid that the concentrations are taken.
    z_phy = np.linspace(-1000, 0, N)
    
    ## Using diatomn as phy abs, scat coefficients
    ab_wat = abscat(lam, 'water')
    ab_phy = abscat(lam, 'Diat') 
    # ab_syn = abscat(lam, 'Syn')
    
    rrs_irr = np.zeros(len(chlor_dat))
    
    for k, chlor in enumerate(chlor_dat):
        ## Constant profile.
        if use_art_phy:
            ## negative
            loc = -70
            width = 40
            phy_prof = OI.artificial_phy_prof(z_phy, loc, width, chlor/2)
        else: 
            phy_prof = np.full(N,chlor)
        
        ## Creating the Phy object.
        phy = OI.Phy(z_phy, phy_prof, ab_phy[0], ab_phy[1])
        ## Running OI
        Eu_surf = OI.ocean_irradiance(z_phy[0], Ed0, Es0, Euh, ab_wat, 
                                      phy = phy, N = 30, 
                                      pt1_perc_zbot = True)[2][-1]
        ## Running R_RS
        rrs_irr[k] = OIR.R_RS(Ed0, Es0, Eu_surf)
        
    return rrs_irr


def Plot_Val_Comparison_Global(longitude, latitude, val1, val2, units, title1, title2):
    
    """
    This plots the given value as a scatter plot of color mapped points 
    on a Cartopy module projection. 

    Parameters
    ----------
    longitude : 1-D Array [N]
        Longitude values. 
    latitude : 1-D Array [N]
        latitude values.
    val1 : 1-D Array [N]
        The first value to be plotted. 
    val2 : 1-D Array [N]
        The second value to be plotted. 
    units : String
        The units of the value being plotted.
    title1 : String
        The title of the first plot, the value being plotted. 
    title2 : String
        The title of the second plot, the value being plotted. 
        
    Returns
    -------
    fig : matplotlib Figure instance
        Figure corresponding to this plot.
    ax1 : matplotlib axes instance
        Axes corresponding to this plot.

    """
    
    
    fig = plt.figure(figsize=(8,6))
    cbar_shrink=.5
    s=15
    vmax = max(val1.max(), val2.max())
    
    ## The first value plot
    ax1 = fig.add_subplot(121, projection=ccrs.Robinson())
    ax1.add_feature(cfeature.COASTLINE)
    ax1.gridlines()
    im1 = ax1.scatter(longitude, latitude, c=val1, s=s, cmap='nipy_spectral', 
                     transform=ccrs.PlateCarree(), vmax=vmax )
    # fig.colorbar(im1, ax=ax1, shrink=cbar_shrink, label = units)
    ax1.set_title(title1)
    
    ## The second value plot
    ax2 = fig.add_subplot(122, projection=ccrs.Robinson())
    ax2.add_feature(cfeature.COASTLINE)
    ax2.gridlines()
    im2 = ax2.scatter(longitude, latitude, c=val2, s=s, cmap='nipy_spectral', 
                     transform=ccrs.PlateCarree(), vmax=vmax )
    fig.colorbar(im1, ax=[ax1,ax2], shrink=cbar_shrink, label = units)
    ax2.set_title(title2)
    
    return fig,ax1,ax2


def Plot_Val_Distibution_Comparison(x,y,xlabel,ylabel,title,lim=None):
    
    fig, ax = plt.subplots()

            
    ax.plot(x, y,'o', fillstyle='none')
    ax.plot(x, x, 'k')
    if lim == None: 
        lim = max(x.max(), y.max())
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    
    
    return 

if __name__ == '__main__':
    file = 'C:/Users/miles/RESEARCH/Current_Work/Ocean_Irradiance/ocean_irradiance_studies/Ocean_Color_Data/1629913643435826_chlor_a.csv'
    
    ## Copied from header of csv file. Might be better to read from file. 
    field_names = "id,latitude,longitude,date_time,cruise,seawifs_filename," +\
    "seawifs_es_error,seawifs_pixel_total,seawifs_tdiff,seawifs_solz,seawifs_senz," +\
    "seawifs_cv,seawifs_windspeed,seawifs_chlor_a,insitu_chlor_a_type,insitu_chlor_a,"+\
    "seawifs_rrs412,seawifs_rrs443,seawifs_rrs490,seawifs_rrs510,seawifs_rrs555,seawifs_rrs670"
    
    ## Reading in csv file
    skip_header = 29
    data_dict = RCD.Make_Data_Dict(file, field_names, skip_header)    

    
    ## Checking out data and changing array type from 'object' to 'float64'
    lat = data_dict['latitude'].astype('float64')
    long = data_dict['longitude'].astype('float64')
    seawifs_chlor = data_dict['seawifs_chlor_a'].astype('float64')
    insitu_chlor = data_dict['insitu_chlor_a'].astype('float64')
    ## 443 nm wavelength rrs
    seawifs_rrs_b_443 = data_dict['seawifs_rrs443'].astype('float64')
    ## 490 nm wavelength rrs
    seawifs_rrs_b_490 = data_dict['seawifs_rrs490'].astype('float64')
    # 555 nm wavelength is taken to be green.
    seawifs_rrs_g = data_dict['seawifs_rrs555'].astype('float64')
    
    
    ## Plotting Data.
    units = 'log(mg chl_a m^-3)'
    # fig,ax1 = Scatter_Plot_Data_Global(long, lat,np.log(seawifs_chlor), units, 'seawifs_chlor_a')
    # Scatter_Plot_Data_Global(long, lat, np.log(insitu_chlor), units, 'insitu_chlor_a')
    
    ## Making the west coast mask. 
    lat_bound = [25,45]
    long_bound = [-135, -105]
    # mask = Mask(long_bound, lat_bound, long, lat, plot=True, ax =ax1)
    mask = Mask(long_bound, lat_bound, long, lat, plot=False)
    
    ## Masking data
    lat_m = lat[mask == 1]
    long_m = long[mask == 1]
    seawifs_chlor_m = seawifs_chlor[mask == 1]
    insitu_chlor_m = insitu_chlor[mask == 1]
    seawifs_rrs_b_443_m = seawifs_rrs_b_443[mask == 1]
    seawifs_rrs_b_490_m = seawifs_rrs_b_490[mask == 1]
    seawifs_rrs_g_m = seawifs_rrs_g[mask == 1]
    
    
    ## Irradiance component of study
    ## ------------------------------
    ## First to try with constant chl_a concentration.
    
    ## Blue and Green wavelengths.

    lam_b_443 = 443
    lam_b_486 = 486
    lam_g = 551
    
    ## The chlorophyll data used for this.
    # chlor_dat = insitu_chlor_m
    chlor_dat = seawifs_chlor_m 
    


    ## RRS Comparisons
    ## ----------------
    rrs_units='sr^-1'
    
    # ##443
    rrs_irr_b_443_m = R_RS_Irr(lam_b_443, chlor_dat, use_art_phy=(False))
    # Plot_Val_Comparison_Global(long_m, lat_m, seawifs_rrs_b_443_m, rrs_irr_b_443_m, rrs_units, 
    #                             'seawifs_rrs_b_443_m', f'rrs_irr_b_m [{lam_b_443} nm]')
    # Plot_Val_Distibution_Comparison(seawifs_rrs_b_443_m, rrs_irr_b_443_m, 'Seawifs RRS',  'Irradiance RRS', "RRS [443nm] ")
    
    ## 486
    rrs_irr_b_486_m = R_RS_Irr(lam_b_486, chlor_dat, use_art_phy=(False))
    # Plot_Val_Comparison_Global(long_m, lat_m, seawifs_rrs_b_490_m, rrs_irr_b_486_m, rrs_units, 
    #                             'seawifs_rrs_b_490_m', f'rrs_irr_b_m [{lam_b_486} nm]')
    # Plot_Val_Distibution_Comparison(seawifs_rrs_b_490_m, rrs_irr_b_486_m, 'Seawifs RRS',  
    #                                 'Irradiance RRS', "RRS [490nm] ")
    
    # ## Green 
    rrs_irr_g_m = R_RS_Irr(lam_g, chlor_dat, use_art_phy=(False))
    # Plot_Val_Comparison_Global(long_m, lat_m, seawifs_rrs_g_m, rrs_irr_g_m, rrs_units, 
    #                             'seawifs_rrs_g_m', f'rrs_irr_g_m [{lam_g} nm]')
    # Plot_Val_Distibution_Comparison(seawifs_rrs_g_m, rrs_irr_g_m, 'Seawifs RRS',  
    #                                 'Irradiance RRS', "RRS Green ")
    
    
    
    
    ## Plotting the rrs irr ratio of blue to green 
    rrs_ratio_443_irr = rrs_irr_b_443_m / rrs_irr_g_m
    rrs_ratio_486_irr = rrs_irr_b_486_m / rrs_irr_g_m
    rrs_ratio_443_seawifs = seawifs_rrs_b_443_m /seawifs_rrs_g_m
    rrs_ratio_490_seawifs = seawifs_rrs_b_490_m /seawifs_rrs_g_m
    
    Plot_Val_Distibution_Comparison(rrs_ratio_443_seawifs, rrs_ratio_443_irr, 
                                    'Seawifs RRS Ratio', 'Irradiance RRS Ratio', 
                                    'RRS Ratio Distribution [Blue = 443 nm]')
    Plot_Val_Distibution_Comparison(rrs_ratio_490_seawifs, rrs_ratio_486_irr, 
                                    'Seawifs RRS Ratio', 'Irradiance RRS Ratio', 
                                    'RRS Ratio Distribution [Blue = 486nm]')

    
    ## chl_a comparisons
    ## -----------------
    
    ## OCx chlorophyll concentrations.
    chl_irr_443 = OIR.OCx_alg(rrs_irr_b_443_m, rrs_irr_g_m)
    # chl_irr_486 = OIR.OCx_alg(rrs_irr_b_486_m, rrs_irr_g_m, lam_b_486, lam_g)
    chl_irr_486 = OIR.OCx_alg(seawifs_rrs_b_443_m, seawifs_rrs_g_m)
    xlabel = 'Chl_a Data [mg chl_a m^-3]'
    ylabel = 'Irradiance Chl_a [mg chl_a m^-3]'
    
    
    Plot_Val_Distibution_Comparison(seawifs_chlor_m, chl_irr_443, 
                                    'Seawifs chl_a', 'Irradiance chl_a', 
                                    'Chl_a comparison [Blue = 443nm]', lim =3)


    



    
    
    
    
    
    
    
    
    
    
    
    
    
    