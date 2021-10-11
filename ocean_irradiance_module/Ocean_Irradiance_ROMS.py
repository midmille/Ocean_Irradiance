# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 07:39:39 2020

@author: Miles Miller
"""
"""
ROMS wrapper for Ocean_Irradiance
"""

import numpy as np
from netCDF4 import Dataset 
## can not import seapy at this point
# import os
## weird thing added because I am missing this value for some reason
# os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share"
# import seapy 
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module import Ocean_Irradiance
import os
import sys
import pickle


def OCx_alg(R_rs_b, R_rs_g) :
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
    
    
    OCx_sol = OCx_alg(R_rs_b, R_rs_g)
    # CI_alg_sol = ocean_color.CI_alg(R_rs_r, R_rs_g, R_rs_b, lam_r, lam_g, lam_b)
    
    return OCx_sol


def Ocean_Irradiance_Field(mask, ab_wat, ab_diat, ab_syn, chl_diatom, chl_nanophyt, 
                  z_r, E_d_0, E_s_0, E_u_h, coefficients, N=30, pt1_perc_zbot = True,
                  method='bvp'):
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
    coeffcients : Tuple, length == 5. 
        Coefficients taken from Dutkiewicz such as the average of cosines and sines. 
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
    ## The size four fourth dimension is for the irradiance field
    ## Ed, Es, Eu, z in that order.
    if method == 'shoot_down' or method == 'shoot_up' or method =='shoot_fp':
        irr_arr = np.zeros((N,nyi,nxi,4)) * np.nan 
    if method == 'Dut':
        irr_arr = np.zeros((N-1,nyi,nxi,4)) * np.nan 
    count = 0
    for j in range(nyi): 
        for i in range(nxi):
            ##land is a zero, only computes it for water
            if mask[j,i] == 1: 
                print("{} out of {}".format(count, (nyi*nxi)))
                count += 1

                ## ROMS vertical grid for this index 
                z_r0 = z_r[:,j,i] 
                # z_w0 = z_w[:,j,i] 
                hbot = z_r0[0]
                assert len(chl_diatom) == len(z_r0)
                
                phy_profs = np.zeros((len(z_r0),2))
                phy_profs[:,0] = chl_diatom[:,j,i]
                phy_profs[:,1] = chl_nanophyt[:,j,i]
                
                a = np.array([ab_diat[0], ab_syn[0]])
                b = np.array([ab_diat[1], ab_syn[1]])
                
                phy = Ocean_Irradiance.Phy(z_r0, phy_profs, a, b)    
                
                if method == 'shoot_down':
                    
                    ocean_irr_sol = Ocean_Irradiance.ocean_irradiance(hbot,E_d_0,E_s_0,E_u_h,
                                                                      ab_wat, coefficients, phy, N=N, 
                                                                      pt1_perc_zbot = pt1_perc_zbot)
                elif method == 'shoot_up':
                    
                    ocean_irr_sol = Ocean_Irradiance.ocean_irradiance_shoot_up(hbot,E_d_0,E_s_0,E_u_h,
                                                                      ab_wat, coefficients, phy, N=N, 
                                                                      pt1_perc_zbot = pt1_perc_zbot)
                elif method == 'shoot_fp':
                    
                    ocean_irr_sol = Ocean_Irradiance.ocean_irradiance_shoot_fp(hbot,E_d_0,E_s_0,E_u_h,
                                                                      ab_wat, coefficients, phy, N=N, 
                                                                      pt1_perc_zbot = pt1_perc_zbot)
                elif method == 'Dut':
                    ocean_irr_sol = Ocean_Irradiance.ocean_irradiance_dutkiewicz(hbot,E_d_0,E_s_0,E_u_h,
                                                                      ab_wat, coefficients, phy, N=N, 
                                                                      pt1_perc_zbot = pt1_perc_zbot)
                ## Ed
                irr_arr[:,j,i,0] = ocean_irr_sol[0]
                ## Es
                irr_arr[:,j,i,1] = ocean_irr_sol[1]
                ## Eu
                irr_arr[:,j,i,2] = ocean_irr_sol[2]
                ## z_irr
                irr_arr[:,j,i,3] = ocean_irr_sol[3]

                
    return irr_arr


def Irradiance_Run(R_nc, PI, nstp, save_dir, save_file, method):
    
    
    ## The name of the file that the python Eu dict will be saved to as pickle.
    # save_file = f'irr_dict_nstp_{nstp}.p'
    # save_dir = 'Irr_Field_Out'
    save_path = f'{save_dir}/{save_file}'
    ## Python calculated Eu at surface 
    ##---------------------------------
    ## Checking if save file exists
    ## if it doesn't exist then redo calculation
    if os.path.exists(save_path) == False:
        ## User input required so not to overwrite or redo unnecessarily... 
        y_or_n = input('File does not exist, continue with calculation? [y/n] ')
        if y_or_n == 'n': 
            sys.exit('Stopping...')
        elif y_or_n == 'y':
            print('ok, statrting calculations...')
            
        #mask = np.ones((R_nc.nyi, R_nc.nxi))
        irr_field_py = {}
        for lam in R_nc.wavelengths:
            print('Current Wavelength:', lam)
            irr_field_py[lam] = Ocean_Irradiance_Field(
                                              R_nc.maskr, 
                                              abscat(lam, 'water'), 
                                              abscat(lam, 'Diat'), 
                                              abscat(lam, 'Syn'), 
                                              R_nc.chl_diatom[nstp,:,:,:], 
                                              R_nc.chl_nanophyt[nstp,:,:,:], 
                                              R_nc.z_r[nstp,:,:,:], 
                                              PI.Ed0, 
                                              PI.Es0, 
                                              PI.Euh,
                                              PI.coefficients,
                                              N= 100, 
                                              method = method)
        
        pickle.dump(irr_field_py, open(save_path, "wb"))
        print('Python calculation complete and saved')
        
    ## if the save file does exist then just load it. gity 
    elif os.path.exists(save_path) == True:
        print(f'Irradiance save file exists! Loading python calculated irradiance field from file "{save_file}"...')
        irr_field_py = pickle.load(open(save_path,'rb'))
        print('Yay, file loaded :)')


# def plot_ocean_color():
#     import matplotlib.pyplot as plt
#     import matplotlib as mpl
    
#     plt.figure()
#     plt.pcolormesh( ocean_color, cmap="nipy_spectral", norm = mpl.colors.LogNorm(), vmin=.01, vmax=66)#, cmap='nipy_spectral')
#     plt.colorbar(label='Chl-a')
#     plt.ylabel('Degrees Lattitude')
#     plt.xlabel('Degress Longitude')
#     plt.title('Ocean Color')
#     plt.show()
    
    
if __name__ == '__main__':
    
    import argparse
    ## User Modules
    from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    parser = argparse.ArgumentParser(description='Ocean Irradiance ROMS Wrapper')
    parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    parser.add_argument('save_file_name', help = "The name of the pickle file in which the result will be stored" )
    parser.add_argument('save_dir', help = "The name and path of the direcotry in which the result will be stored" )
    parser.add_argument('method', help = "Methods are: Dut, shoot_down, shoot_up, shoot_fp." )
    # parser.add_argument('dest_file', help='Path to Destination Directory. Saved as pickle')
    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    args = parser.parse_args()
    
    print(args.file)

    PI = Param_Init()
    R_nc = ROMS_netcdf(args.file)

    ## updating default to wavelengths necessary for the OCx_algorithim
    R_nc.wavelengths = [443,551]

    ## settting time step index
    time_step_index = 1
    
    Irradiance_Run(R_nc, PI, time_step_index, args.save_dir, args.save_file_name, args.method)
    
    
    
    
    
    
    ## THIS IS DEPRICATED!!! JUST ADD TO R_nc OBject.... R_nc does not inclued NetCDF object any more 
    ## This is due to difficulty with pickling.
    # lon_rho = R_nc.roms_nc.variables['lon_rho'][:]
    # lat_rho = R_nc.roms_nc.variables['lat_rho'][:]

    # ## Calculating Eu at Surface dictionary. 
    # #--------------------------------------------------------------------------
    # Eu_surface_dict = {}
    # for lam in PI.wavelengths: 
        
    #     ##eyeballed absorbtion and scattering from Dut. 2015
    #     ab_wat = abscat(lam, 'water')
    #     ab_diat = abscat(lam, 'Diat') 
    #     ab_syn = abscat(lam, 'Syn')
        
    #     Eu_surface_dict[lam] = Ocean_Irradiance_Field(lam, R_nc.maskr, ab_wat, ab_diat, ab_syn, 
    #                                          R_nc.diatom, R_nc.nanophyt, R_nc.z_w, 
    #                                          R_nc.z_r, PI.Ed0, PI.Es0, PI.Euh,
    #                                          PI.coefficients, PI.Chl2NL, PI.Chl2NS)[-1, :,:, 2]
    
    # ## Ocean Color Calculation
    # #--------------------------------------------------------------------------
    # ocean_color = ocean_color_sol(Eu_surface_dict, PI.Ed0, PI.Es0)
    # pickle.dump(ocean_color, open("ocean_color.p", "wb"))
    
    # ## Plotting
    # #--------------------------------------------------------------------------
    # if args.plot:

    #     plot_ocean_color()
    
