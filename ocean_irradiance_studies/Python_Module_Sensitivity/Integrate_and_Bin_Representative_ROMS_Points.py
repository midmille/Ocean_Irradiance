# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 06:25:01 2021

@author: Miles Miller

--> This script helps to search for a limited set of points that reflect the different 
phytoplankton concentration profiles within the ROMS output domain.
 
--> It functions basically by numerically integrating vertically over all profiles 
within the domainand then placing the resulting integral solutions into bins 
and finding the index of one random point from within each bin. 

--> This helps limit the field of calculation in the case of computationally 
costly sensitivity studies and lessen time for computation. 

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
import PARAMS

## Importing other modules. 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import numpy as np
import random 


def Rep_ROMS_Phy_Prof_Loc(nxi, nyi, chl_diatom, chl_nanophyt,
                          z_r, n_bins, db =20, plot_flag=True) :
    """
    This function calculates the index locations of a set of points representative 
    of different phytoplankton profiles. 

    Parameters
    ----------
    nxi : Float
        The number of points in the x direction, number of columns. 
    nyi : Float
        The number of points in the y direction, number of rows. 
    chl_diatom : 3-D Array
        The diatom chl concentration array. 
    chl_nanophyt : 3-D Array
        The diatom chl concentration array.
    z_r : 3-D Array
        The rho grid. The cell centers. 
    n_bins : Float
        The number of bins to seperate into. Also the number of representative 
        indexes to be returned. 
    db : Float, optional
        The boundary distance in indices. The default is 20.
    plot_flag : Boolean, optional
        To plot or not. Plots the position of representative points on a contour of 
        the integrals. The default is True.

    Returns
    -------
    bin_point_index: 2-D Array
        The index positions of the representative points on the unbounded ROMS 
        domain. 
    

    """
    
    
    def Integrate_Phy_Prof(nyi_db, nxi_db, phyt_profiles, zarr_phy) :
        """
        Integrates over the vertical phytoplankton profiles within the horizontal domain. 
        Each point in the horizontal ROMS domain is then represented by a
        float of the integral. 
        

        Parameters
        ----------
        phyt_profiles : 3-D Array
            The phytoplankton concentration profiles. [z,y,x]
        zarr_phy : 3-D Array
            The corresponding coordinates of the phytoplankton concentrations. 

        Returns
        -------
        integral_phy : 2-D Array
            The vertical integral of the phytoplankton profiles for each point 
            in the horizontal domain. 

        """
        ## Empty integral array of same shape as the phyt_profile array.
        integral_phy = np.zeros_like(phyt_profiles[0,:,:])
        
        ## Looping over the horizontal domain. 
        for k in range(nyi_db) :
            for j in range(nxi_db) :
                ## Integrating over the vertical direction of array. 
                integral_phy[k,j] = np.trapz(phyt_profiles[:,k,j], zarr_phy[:,k,j]) 
        
        return integral_phy
    
    
    def Find_Random_Point_In_Bin(integral_map, bin_edges, points_in_bins) :
        """
        This finds the index of a random point that lies within each integral bin. 

        Parameters
        ----------
        integral_map : 2-D Array
            The horizontal map of integral values within bounded domain. 
        bin_edges : 1-D Array
            The value of the bin edges. Should be of length one greater than the 
            number of bins. 
        points_in_bins : 1-D Array
            The number of points in each bin. The values of each bar in the histogram.

        Returns
        -------
        bin_point_index : 2-D Array
            The index of the representative points on the bounded domain. 

        """
        ## The index location of the chosen points within the bounded ROMS domain.
        ## The 2 in the second dimension is for (y,x)
        bin_point_index = np.zeros((len(points_in_bins), 2))
        ## Randomizing the indexes so that the points are taken in different locations in domain. 
        rand_nyi_arr = np.arange(0, nyi_db)
        rand_nxi_arr = np.arange(0, nxi_db)
        random.shuffle(rand_nyi_arr)
        random.shuffle(rand_nxi_arr)
        
        ## Looping over the bins.
        for k in range(len(points_in_bins)) : 
            ## If there are no points in the bin then return the location as NaN.
            if points_in_bins[k] == 0 : 
                bin_point_index[k,:] = np.NaN
            ## If there are point sin the bin then loop over the domain
            ## to find the index of a point that is whithin that bin. 
            if points_in_bins[k] > 0 :
                ## Loop over the horizontal domain.
                for yi in rand_nyi_arr :
                    for xi in rand_nxi_arr : 
                        ## If the value of the integral lies within the current bin 
                        ## set the bin_point_index for that integral point to be 
                        ## the current index and then break out of both loops. 
                        if bin_edges[k] < integral_map[yi,xi] < bin_edges[k+1]:
                            bin_point_index[k,:] = [yi,xi]
                            break 
                    else: 
                        continue 
                    break 
                
        
        return bin_point_index

    ## Taking only points in arays within the boundary of db. 
    chl_diatom_db = chl_diatom[:,db:nyi-db, db:nxi-db]
    chl_nanophyt_db = chl_nanophyt[:, db:nyi-db, db:nxi-db]
    
    ## The number of points in y,x directions in boundary domain
    nyi_db = len(chl_diatom_db[0,:,0])  
    nxi_db = len(chl_diatom_db[0,0,:])

    ## The z grid at which the phytoplankton concentrations are located.
    ## Spliced to the correct domain within the boundary db. 
    zarr_phy = np.copy(z_r[:,db:nyi-db,db:nxi-db])  
    
    ## Integrating over vertical profiles in bounded domain as a metric of the
    ## concentration profile structure.
    diatom_integrate = Integrate_Phy_Prof(nyi_db, nxi_db, chl_diatom_db, zarr_phy)
    nanophyt_integrate = Integrate_Phy_Prof(nyi_db, nxi_db, chl_nanophyt_db, zarr_phy)
    ##Total integral of points
    integral_map = diatom_integrate + nanophyt_integrate
    ## Log base 10 of integral field
    log10_integral_map = np.log10(integral_map)
    ## Binning distribution 
    points_in_bins, bin_edges = np.histogram(log10_integral_map, n_bins)
    
    ## Getting the y,x index of random points that reside within each bin. 
    bin_point_index = Find_Random_Point_In_Bin(log10_integral_map, bin_edges, points_in_bins)
    ## Adding the boundary value back so that the indexes represent the indexes on 
    ## the unbounded ROMS grid such that and index of (0,0) on the bounded domain
    ## would become (db,db) on the unbounded.
    bin_point_index = bin_point_index + db 
    
    if plot_flag == True : 
        ## Integrating over vertical profiles in FULL non-bounded domain as a metric of the
        ## concentration profile structure.
        diatom_integrate_full = Integrate_Phy_Prof(nyi, nxi, chl_diatom, z_r)
        nanophyt_integrate_full = Integrate_Phy_Prof(nyi, nxi, chl_nanophyt, z_r)
        ##Total integral of points in full domain
        integral_map_full = diatom_integrate_full + nanophyt_integrate_full
        ## Log base 10 of integral field
        log10_integral_map_full = np.log10(integral_map_full)
        ##plotting diatom 
        fig, ax = plt.subplots()
        ax.pcolormesh(log10_integral_map_full)
        # for k in range(len(bin_point_index)) :
        ax.plot(bin_point_index[:,1], bin_point_index[:,0], 'ro', markersize = 4, label='Location of Chosen Profiles')
        ax.set_ylabel('nyi')
        ax.set_xlabel('nxi')
        ax.set_title(r'$\log _{10} (\int \mathrm{diatom} + \int \mathrm{nano} )$')
        ax.legend(title=f'N Profiles = {n_bins}')
    
    return bin_point_index
  
if __name__ == '__main__':
    
    ## Creating a parameter object. 
    file = '/Users/Miles Miller/da_fwd_002.nc'
    p = PARAMS.Param_Init(file)
    ## The number of representative points.
    n_bins = 20
    Rep_ROMS_Phy_Prof_Loc(p.nxi, p.nyi, p.chl_diatom, p.chl_nanophyt,
                          p.z_r, n_bins, db =20, plot_flag=True)













