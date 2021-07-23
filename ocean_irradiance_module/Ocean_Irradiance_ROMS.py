# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 07:39:39 2020

@author: Miles Miller
"""
"""
ROMS wrapper for Ocean_Irradiance
"""

import numpy as np
import Ocean_Irradiance
from netCDF4 import Dataset 
import os
## weird thing added because I am missing this value for some reason
os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share"
import seapy 
from absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
import pickle

def make_roms_grid_from_netcdf(file):
    """
    Retreives the ROMS grid from file. 
    

    Parameters
    ----------
    file : string 
        The netcdf file output from roms. 

    Returns
    -------
    z_r : Array
        The rho grid. The cell centers 
    z_w : Array 
        The w grid. The cell edges. 

    """
        
    roms_data = Dataset(file, 'r')
    
    h = roms_data.variables['h'][:]
    hc = roms_data.variables['hc'][:]
    s_rho = roms_data.variables['s_rho'][:]
    s_w = roms_data.variables['s_w'][:]
    Cs_r = roms_data.variables['Cs_r'][:]
    Cs_w = roms_data.variables['Cs_w'][:]
    zeta = roms_data.variables['zeta'][0,:]
    vtransform = roms_data.variables['Vtransform'][:]
    
    z_r = seapy.roms.depth(vtransform, h, hc, s_rho, Cs_r,zeta)[:,:,:] ##the rho grid 
    z_w = seapy.roms.depth(vtransform, h, hc, s_w, Cs_w,zeta)[:,:,:]  ##the w-grid at eddges 
    
    return z_r, z_w


def OCx_alg(R_rs_b, R_rs_g, lam_b, lam_g) :
    """
    

    Parameters
    ----------
    R_rs_b : Float
        R_rs in the blue wavelength
    R_rs_g : Float 
        R_rs in the green wavelength
    lam_b : Float
        Wavelength of blue
    lam_g : Float 
        Wavelength of green 

    Returns
    -------
    chl_a : Float 
        Value of chl_a

    """
    ## the a values for OCx 0.2228	-2.4683	1.5867	-0.4275	-0.7768
    a0 = .2228
    a1 = -2.4683
    a2 = 1.5867
    a3 = -0.4275
    a4 = -.7768
    
    
    log_10_chl_a = a0 ##log base 10 of chl-a conc. 
    for a,i in zip([a1,a2,a3,a4],[1,2,3,4]): 
        
        log_10_chl_a = log_10_chl_a + a * ((np.log10((R_rs_b)/(R_rs_g)))**i)
    chl_a = 10**(log_10_chl_a)
    return chl_a 


def R_RS(E_d_0, E_s_0, E_u_surface ) : 
    """
    Calculates the remotely sensed reflectance following Dutkiewicz (2015) equations 
    5 and 6. 

    Parameters
    ----------
    E_d_0 : Float 
        Initial value of downward direct irradiance. 
    E_s_0 : Float 
        Initial value of downward diffuse irradiance. 
    E_u_surface : Float
        Solution for upwelling irradiance at the surface. 

    Returns
    -------
    R_rs : Float
        The remotely sensed reflectance. 

    """
    
    R = (E_u_surface) / (E_d_0 + E_s_0) ## This is Equation 5 
    
    ##We are assuming that Q =about 4 (Dut. 2015 Eqn. 6)
    Q = 4 
    R_rs = R / Q
    
    return R_rs 


def ocean_color_sol(Eu_sol_dict, E_d_0, E_s_0): 
    """
    Calculates the ocean color solution for the surface value of chl_a

    Parameters
    ----------
    Eu_sol_dict : Dict, [key=wavelength, value=Eu_array]
        A dictionary with the solution for wavelength dependent upwelling irradiance 
    E_d_0 : Float
        Initial value of the downward direct irradiance. 
    E_s_0 : Float
        Initial value of the downward diffuse irradiance. 

    Returns
    -------
    OCx_sol : Float
        The solution for chl_a concentration at the surface. 

    """
    ##the closest our wavelenghts are to the values of ocean color lam
    ## 443(ours) = 443(theirs); 551(ours) = 555(theirs); 671(ours) = 670(theirs) 
    lam_b = 443
    lam_g = 551
    # lam_r = 671
    
    ##getting the R_rs from upwelling and the initial conditiongs
    
    R_rs_b = R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_b]) ##R_rs as a function of the blue wavelength 
    R_rs_g = R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_g]) ##green 
    # R_rs_r = R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_r]) ##red 
    
    
    OCx_sol = OCx_alg(R_rs_b, R_rs_g, lam_b, lam_g)
    # CI_alg_sol = ocean_color.CI_alg(R_rs_r, R_rs_g, R_rs_b, lam_r, lam_g, lam_b)
    
    return OCx_sol


def Eu_at_surface(mask, ab_wat, ab_diat, ab_syn, chl_diatom, chl_nanophyt, 
                  z_r, E_d_0, E_s_0, E_u_h, N=30, pt1_perc_zbot = True):
    """
    

    Parameters
    ----------
    mask : 2-D Array
        Land is taken to be zero, water is one. 
    ab_wat : Tuple, (a,b)
        The absorbtion and scattering coefficients for water. 
    ab_diat : Tuple, (a,b)
        The absorbtion and scattering coefficients for diatoms
    ab_syn : Tuple, (a,b)
        The absorbtion and scattering coefficients for syn. 
    chl_diatom : 3-D Array
        Diatom concentrations from ROMS in chlorophyll.
    chl_nanophyt : 3-D Array
        Nanophytoplankton concentrations from ROMS in chlorophyll.
    z_r : 3-D Array
        The ROMS grid rho point located at cell centers. 
    E_d_0 : Float
        Initial value for downward direct irradiance. 
    E_s_0 : Float
        Initial value for downward diffuse irradiance. 
    E_u_h : FLoat
        Bottom boundary condition on upwelling irradiance. 
    N : Float, default is 30
        The number of layers in the logarithmic grid. 
    pt1_perc_zbot : Boolean, default is True
        True refers to using the .1% light level as the zbot so long as that the magnitude 
        of the .1% light level is smaller than the magnitude of hbot. False refers
        to just using given hbot as zbot. 

    Returns
    -------
    Eu_arr : 2-D Array
        The array of surface values of upwelling irradiance. 

        

    """

    nyi,nxi = np.shape(mask)
    
    ##Not a Number array so the land is Nan, not 0, helps make land white in pcolormesh
    Eu_arr = np.zeros((nyi,nxi)) * np.nan 
    count = 0
    for j in range(nyi): 
        for i in range(nxi):
            ##land is a zero, only computes it for water
            if mask[j,i] == 1: 
                print("{} out of {}".format(count, (nyi*nxi)))
                count += 1
                ##from Nitogen to Chl
                chl_diatom = (diatom[:,j,i]*Chl2NL)
                chl_nanophyt = (nanophyt[:,j,i]*Chl2NS)
                
                ## ROMS vertical grid for this index 
                z_r0 = z_r[:,j,i] 
                # z_w0 = z_w[:,j,i] 
                zbot = z_r0[0]
                assert len(chl_diatom) == len(z_r0)
                
                phy_profs = np.zeros((len(z_r0),2))
                phy_profs[:,0] = chl_diatom
                phy_profs[:,1] = chl_nanophyt
                
                a = np.array([ab_diat[0], ab_syn[0]])
                b = np.array([ab_diat[1], ab_syn[1]])
                
                phy = Ocean_Irradiance.Phy(z_r0, phy_profs, a, b)    
                
                ocean_irr_sol = Ocean_Irradiance.ocean_irradiance(zbot,E_d_0,E_s_0,E_u_h,
                                                                ab_wat,phy, N=N, 
                                                                pt1_perc_zbot = pt1_perc_zbot)
                
                Eu_arr[j,i] = ocean_irr_sol[2][-1]

                
    return Eu_arr


if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Ocean Irradiance ROMS Wrapper')
    parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    # parser.add_argument('dest_file', help='Path to Destination Directory. Saved as pickle')
    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    args = parser.parse_args()
    
    file = args.file
    
    ###################################PARAMS##################################
    # file = '/Users/Miles Miller/da_fwd_002.nc'

    Ed0 = .7
    Es0 = 1 - Ed0
    Euh = 0 
    
    # wavelengths = [410,443,486,551,638,671] 
    ## wavelengths necessary for the OCx_algorithim
    wavelengths = [443,551]
    
    Chl2NL = 1.59
    Chl2NS = .7950
    
    time_step_index = 0 
    ###################################PARAMS END##############################
    
    ## ROMS
    #-------------------------------------------------------------------------
    roms_nc = Dataset(file, 'r')
    roms_date = seapy.roms.num2date(nc = roms_nc)[0]
    
    lon_rho = roms_nc.variables['lon_rho'][:]
    lat_rho = roms_nc.variables['lat_rho'][:]
     
    z_r,z_w = make_roms_grid_from_netcdf(file) ##get the roms grid from the roms_grid module
    
    ocean_time = roms_nc.variables['ocean_time'][:] ##times in seconds since 1999, 01/01, 00:00:003
    
    nti,nzi,nyi,nxi = roms_nc.variables["diatom"].shape ## shaped of the diatom output 
    
    maskr = roms_nc.variables['mask_rho'][:] ##land or not 
    
    roms_date = seapy.roms.num2date(nc = roms_nc)
    
    diatom = roms_nc.variables['diatom'][time_step_index,:,:,:] 
    nanophyt = roms_nc.variables['nanophytoplankton'][time_step_index,:,:,:]
    
    nti,nzi,nyi,nxi = roms_nc.variables["diatom"].shape ##shaped of the diatom output

    ## Calculating Eu at Surface dictionary. 
    #--------------------------------------------------------------------------
    Eu_surface_dict = {}
    for lam in wavelengths: 
        
        ##eyeballed absorbtion and scattering from Dut. 2015
        ab_wat = abscat(lam, 'water')
        ab_diat = abscat(lam, 'Diat') 
        ab_syn = abscat(lam, 'Syn')
        
        Eu_surface_dict[lam] = Eu_at_surface(lam, maskr, ab_wat, ab_diat, ab_syn, 
                                             diatom, nanophyt, z_w, z_r, Ed0, Es0, Euh,
                                             Chl2NL, Chl2NS)
    
    ## Ocean Color Calculation
    #--------------------------------------------------------------------------
    ocean_color = ocean_color_sol(Eu_surface_dict, Ed0, Es0)
    pickle.dump(ocean_color, open("ocean_color.p", "wb"))
    
    ## Plotting
    #--------------------------------------------------------------------------
    if args.plot:
        def plot_ocean_color():
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            
            plt.figure()
            plt.pcolormesh(lon_rho, lat_rho, ocean_color, cmap="nipy_spectral", norm = mpl.colors.LogNorm(), vmin=.01, vmax=66)#, cmap='nipy_spectral')
            plt.colorbar(label='Chl-a')
            plt.ylabel('Degrees Lattitude')
            plt.xlabel('Degress Longitude')
            plt.title('Ocean Color')
            plt.show()
        plot_ocean_color()
    
