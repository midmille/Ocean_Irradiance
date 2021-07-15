# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 06:25:01 2021

@author: Miles Miller
"""

## Importing Ocean_Irradiance modules. 
import os
cwd = os.getcwd()
os.chdir('../ocean_irradiance_module')
import Ocean_Irradiance 
import Ocean_Irradiance_ROMS
from absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
os.chdir(cwd)
import wavelength_to_rgb
## Importing other modules. 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import numpy as np
import random 


###############################Representative ROMS points#####################
#----------------------------------------------------------------------------#


def Rep_ROMS_phy_prof_loc(roms_data, time_step_index, nxi, nyi, chl_diatom, chl_nanophyt,
                          z_r, bins, plot_flag=True) :
    
    
    def integrate(phyt_profiles, zarrs, bins=20) :
        
        integrals = np.zeros_like(phyt_profiles[0,:,:])
        
        for k in range(nyi_db) :
            for j in range(nxi_db) :
                integrals[k,j] = np.trapz(phyt_profiles[:,k,j], zarrs[:,k,j]) 
        
        return integrals
    
    def Find_Point_In_Bin(integral_map, bin_edges, points_in_bins) :
        bin_point_index = np.zeros((len(points_in_bins), 2)) #2 for two location indexes 
        
        rand_nyi = np.arange(0, nyi_db)#np.arange(dist_from_boundary,nyi-dist_from_boundary) ##range excludes the boundaries of domain with the 5
        rand_nxi = np.arange(0, nxi_db)#np.arange(dist_from_boundary,nxi-dist_from_boundary)
        random.shuffle(rand_nyi)
        random.shuffle(rand_nxi)
        
        for k in range(len(points_in_bins)) : 
            if points_in_bins[k] == 0 : 
                bin_point_index[k,:] = np.NaN
            if points_in_bins[k] > 0 :
                for yi in rand_nyi :
                    for xi in rand_nxi : 
                        if bin_edges[k] < integral_map[yi,xi] < bin_edges[k+1]:
                            bin_point_index[k,:] = [yi,xi]
                        # else: 
                        #     print('mayday')
        
        return bin_point_index
    
    ## The distance in indexes from boundary. 
    db = 20 
    ## Taking only points in arays within the boundary of db. 
    chl_diatom_db = chl_diatom[:,db:nyi-db, db:nxi-db]
    chl_nanophyt_db = chl_nanophyt[:, db:nyi-db, db:nxi-db]
    
    ## The number of points in y,x directions in boundary domain
    nyi_db = len(chl_diatom_db[0,:,0])  
    nxi_db = len(chl_diatom_db[0,0,:])

    ## The z grid at which the phytoplankton concentrations are located.
    ## Spliced to the correct domain within the boundary db. 
    zarr = np.copy(z_r[:,db:nyi-db,db:nxi-db])  
    
    
    ## Integrating over vertical profiles in domain as a metric of the concentration...
    ## profile structure.
    diatom_integrate = integrate(chl_diatom_db, zarr)
    nanophyt_integrate = integrate(chl_nanophyt_db, zarr)
    ##total integral of points
    integral_map = diatom_integrate + nanophyt_integrate
    ##log base 10 of integral field
    log10_integral_map = np.log10(integral_map)
    ##binning distribution 
    points_in_bins, bin_edges = np.histogram(log10_integral_map, bins)
    
    bin_point_index = Find_Point_In_Bin(log10_integral_map, bin_edges, points_in_bins)
    
    if plot_flag == True : 
        ##plotting diatom 
        fig, ax = plt.subplots()
        ax.pcolormesh(log10_integral_map)
        # for k in range(len(bin_point_index)) :
        ax.plot(bin_point_index[:,1], bin_point_index[:,0], 'ro', markersize = 6, label='Location of Chosen Profiles')
        ax.legend()
    
    return bin_point_index, integral_map, bin_edges, points_in_bins
  
Rep_ROMS_phy_prof_loc(bins=10, plot_flag=True)  