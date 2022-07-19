# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 13:38:54 2020

@author: Miles Miller
"""

"""
ROMS Results and Plot

File: da_fwd_002.nc



"""
from netCDF4 import Dataset 

import matplotlib.pyplot as plt 
import matplotlib as mpl
# import ocean_irradiance_roms_current_edition as irr
import ocean_irradicance_current_edition_2_BC_Eubot_is_zero as irr
# import ocean_irradicance_current_edition_2_BC_Eubot_is_zero_2_edited_for_debug as irr
# import ocean_irradiance_cell_center as irr
from absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
import numpy as np
import wavelength_to_rgb
import ROMS_grid
# import os
# os.environ['PROJ_LIB'] = 'C:/Users/Miles Miller/Anaconda3/envs/pybasemap36/lib/site-packages/mpl_toolkits/basemap' ##not really sure what this does 
import seapy 
import time
import multiprocessing
import nasa_ocean_color_algorithim as ocean_color
import Three_Stream_numerical_derivative as irr_numerical
import ocean_irradiance_bvp as irr_bvp
import ocean_irradiance_numerical_optimization as irr_num_opt

#####################################INITIAL PARAMS###########################
##############################################################################
file = '/Users/Miles Miller/da_fwd_002.nc'
E_d_0 = .7
E_s_0 = 1 - E_d_0

wavelengths = [410,443,486,551,638,671] ##Taken from John Wilkin email regarding sattelite waveelengths
# wavelengths = [410]
roms_data = Dataset(file, 'r')

Chl2NL = 1.59
Chl2NS = .7950

time_step_index = 0 

##NOTE: this takes a little while to compute, make a cell to use other parts of this file 
z_r,z_w = ROMS_grid.make_roms_grid_from_netcdf(file) ##get the roms grid from the roms_grid module

ocean_time = roms_data.variables['ocean_time'][:] ##times in seconds since 1999, 01/01, 00:00:003

nti,nzi,nyi,nxi = roms_data.variables["diatom"].shape ##shaped of the diatom output 

maskr = roms_data.variables['mask_rho'][:] ##land or not 

roms_date = seapy.roms.num2date(nc = roms_data)

zarr = np.linspace(-200,0,201) ##1m res to 200m 

##############################################################################
##############################################################################

#%%

##################################Solving for Eu at surface for grid point####
##############################################################################

def Eu_at_surface(wavelength):
    """
    

    Parameters
    ----------
    wavelength : Integer
        Any wavelength from wavelengths list. 

    Returns
    -------
    Eu_arr : 2-D Array 
        The surface value of the upwelling irradiance. 

    """
    lam = wavelength
    ab_wat = abscat(lam,'water') ##abscat is the file of all the different 
    ab_diat = abscat(lam, 'Diat') ##eyeballed absorbtion and scattering from Dut. 2015
    ab_syn = abscat(lam, 'Syn')
    bad_points_lt = 0
    bad_points_gt = 0
    good_points = 0
    Eu_arr = np.zeros((nyi,nxi)) * np.nan ##Not a Number array so the land is Nan, not 0, helps make land white in pcolormesh
    for j in range(nyi): 
        for i in range(nxi):
            if maskr[j,i] == 1: ##land is a zero, only computes it for water

                diatom = roms_data.variables['diatom'][time_step_index,:,j,i] 
                nanophyt = roms_data.variables['nanophytoplankton'][0,:,j,i] 
                ##Checking if any of the phytoplankton profile is negative 
                ##from Nitogen to Chl
                chl_diatom = (diatom*Chl2NL)  ## Miles Will confirm this.  
                chl_nanophyt = (nanophyt*Chl2NL)
                z_r0 = z_r[:,j,i] ##the verticle array of z_r, z_w
                z_w0 = z_w[:,j,i]
                chl_w_diatom = np.interp(z_w0,z_r0,chl_diatom) ##interpolating the chl to the z_w grid
                chl_w_nanophyt = np.interp(z_w0,z_r0,chl_nanophyt)
                phy_prof_and_ab_tuple_list = [(chl_w_diatom,ab_diat[0], ab_diat[1]), 
                                              (chl_w_nanophyt, ab_syn[0], ab_syn[1])]
                #[(0*z_w0,0,0)]
                three_stream = irr.ocean_irradiance(z_w0, z_w0, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)
                E_limit_upper = 1 ##irradiance limit
                E_limit_lower = 0
                ##irradiance max and min value flag 
                if np.any(three_stream[0] > E_limit_upper) or np.any(three_stream[1] > E_limit_upper) or np.any(three_stream[2] > E_limit_upper):
                    # print("BAD POINT ERROR: Eu, Es, or Ed greater than 1 \n",'wavelength:',lam,'\n index = ', j,',', i )
                    bad_points_gt += 1
                elif np.any(three_stream[0] < E_limit_lower) or np.any(three_stream[1] < E_limit_lower) or np.any(three_stream[2] < E_limit_lower):
                    # print("BAD POINT ERROR: Eu, Es, or Ed less than 0 \n",'wavelength:',lam,'\n index = ', j,',', i )
                    bad_points_lt += 1
                else: 
                    print('good point:',wavelength,':',j, i)
                    good_points += 1 
                # print(wavelength,':',j, i)
                Eu_arr[j,i] = three_stream[-1][-1]
                
    return Eu_arr, bad_points_gt, bad_points_lt, good_points


def Eu_at_surface_zarr_interp(wavelength):
    """
    Calculates the Eu_at_surface with a 1m uniform grid. 

    Parameters
    ----------
    wavelength : Integer
        Any wavelength from wavelengths list. 

    Returns
    -------
    Eu_arr : 2-D Array 
        The surface value of the upwelling irradiance. 

    """
    ab_wat = abscat(wavelength,'water') ##abscat is the file of all the different 
    ab_diat = abscat(wavelength, 'Diat') ##eyeballed absorbtion and scattering from Dut. 2015
    ab_syn = abscat(wavelength, 'Syn')
    Eu_arr = np.zeros((nyi,nxi)) * np.nan ##Not a Number array so the land is Nan, not 0, helps make land white in pcolormesh
    for j in range(nyi): 
        for i in range(nxi):
            if maskr[j,i] == 1: ##land is a zero, only computes it for water
                diatom = roms_data.variables['diatom'][0,:,j,i] 
                nanophyt = roms_data.variables['nanophytoplankton'][0,:,j,i] 
                if np.any(diatom < 0): 
                    print('negative diat concentation')
                neg_phy_index = np.where(np.any(diatom < 0))[0]
                for i in neg_phy_index: 
                    diatom[i] = 0 
                ##from Nitogen to Chl
                chl_diatom = (diatom*Chl2NL)  ## Miles Will confirm this.  
                chl_nanophyt = (nanophyt*Chl2NL)
                z_r0 = z_r[:,j,i] ##the verticle array of z_r, z_w
                if z_w[0,j,i] > zarr[0]:
                    z_w0 = np.linspace(z_w[0,j,i],z_w[-1,j,i], int(z_w[-1,j,i] - z_w[0,j,i]))
                    # z_w0 = np.linspace(z_w[0,j,i],z_w[-1,j,i], int(len(zarr)))
                elif z_w[0,j,i] < zarr[0]:
                    z_w0 = zarr
                chl_w_diatom = np.interp(z_w0,z_r0,chl_diatom) ##interpolating the chl to the z_w grid
                chl_w_nanophyt = np.interp(z_w0,z_r0,chl_nanophyt)
                phy_prof_and_ab_tuple_list = [(chl_w_diatom,ab_diat[0], ab_diat[1]), 
                                              (chl_w_nanophyt, ab_syn[0], ab_syn[1])]
                
                Eu_arr[j,i] = irr.ocean_irradiance(z_w0, z_w0, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)[-1][-1]
                  
    return Eu_arr


def Eu_at_surface_bvp_sol(wavelength) : 
    ab_wat = abscat(wavelength,'water') ##abscat is the file of all the different 
    ab_diat = abscat(wavelength, 'Diat') ##eyeballed absorbtion and scattering from Dut. 2015
    ab_syn = abscat(wavelength, 'Syn')
    Eu_arr = np.zeros((nyi,nxi)) * np.nan ##Not a Number array so the land is Nan, not 0, helps make land white in pcolormesh
    diatom = roms_data.variables['diatom'][0,:,:,:] 
    nanophyt = roms_data.variables['nanophytoplankton'][0,:,:,:] 
    for j in range(nyi): 
        for i in range(nxi):
            if maskr[j,i] == 1: ##land is a zero, only computes it for water
                diatom0 = diatom[:,j,i]
                nanophyt0 = nanophyt[:,j,i]

                ##from Nitogen to Chl
                chl_diatom = (diatom0*Chl2NL)  ## Miles Will confirm this.  
                chl_nanophyt = (nanophyt0*Chl2NL)
                z_r0 = z_r[:,j,i] ##the verticle array of z_r, z_w
                print(wavelength,':',j, i)
                z_w0 = z_w[:,j,i]
                chl_w_diatom = np.interp(z_w0,z_r0,chl_diatom) ##interpolating the chl to the z_w grid
                chl_w_nanophyt = np.interp(z_w0,z_r0,chl_nanophyt)
                phy_prof_and_ab_tuple_list = [(chl_w_diatom,ab_diat[0], ab_diat[1]), 
                                              (chl_w_nanophyt, ab_syn[0], ab_syn[1])]
                
                Eu_arr[j,i] = irr_bvp.ocean_irradiance(z_w0, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)[-1][-1]
                  
    return Eu_arr
    
##Creating a dictionary with the wavelengths as keys, Eu arrays as values. 
##starting the timer 
# start = time.time()
def make_Eu_dict() :
    Eu_sol_dict = {} ##dictionary for the different Eu solutions at each wavelength
    bad_points_gt_dict = {}
    bad_points_lt_dict = {}
    good_points_dict = {}
    for lam in wavelengths:
        Eu_at_surface_sol = Eu_at_surface(lam)
        Eu_sol_dict[lam] = Eu_at_surface_sol[0]
        bad_points_gt_dict[lam] = Eu_at_surface_sol[1]
        bad_points_lt_dict[lam] = Eu_at_surface_sol[2]
        good_points_dict = Eu_at_surface_sol[3]
    
    return Eu_sol_dict, bad_points_gt_dict, bad_points_lt_dict, good_points_dict
    
# if __name__ == '__main__':
#     pool = multiprocessing.Pool(processes = 3)
#     Eu_sol = (pool.map(Eu_at_surface_zarr_interp, wavelengths))
#%%
# Eu_sol_dict = {}
# for k in range(len(Eu_sol)) :
    
#     Eu_sol_dict[wavelengths[k]] = Eu_sol[k]


#%%

################################Wavelength Composite##########################
##############################################################################

def wavelength_composite(Eu_sol_dict):
    """
    

    Parameters
    ----------
    Eu_sol_dict : Dict
        The dictionary of each wavelength to the corresponding solution of 
        Eu from Eu_at_surface function.

    Returns
    -------
    rgb_grid : 3-D array 
        The red, green, and blue value as a float less than one for each point 
        in the roms grid. 

    """
    ##wavelength to rgb  
    lam_rgb_dict = {} 
    for lam in wavelengths: 
        lam_rgb_dict[lam] = np.array(wavelength_to_rgb.wavelength_to_rgb(lam))
    
    r_I_sum = 0
    g_I_sum = 0
    b_I_sum = 0
    I_sum = 0 
    # for j in range(nyi): 
    #         for i in range(nxi):
    for lam in wavelengths: 
        r_I_sum = r_I_sum + Eu_sol_dict[lam] * lam_rgb_dict[lam][0]
        g_I_sum = g_I_sum + Eu_sol_dict[lam] * lam_rgb_dict[lam][1]
        b_I_sum = b_I_sum + Eu_sol_dict[lam] * lam_rgb_dict[lam][2]
        I_sum = I_sum + Eu_sol_dict[lam]
    r_norm = (r_I_sum / I_sum) / 255
    g_norm = (g_I_sum / I_sum) / 255
    b_norm = (b_I_sum / I_sum) / 255
    
    rgb_grid = np.zeros((nyi,nxi,3)) ##three deep for the rgb values
    for j in range(nyi):
        for i in range(nxi):
            rgb_grid[j,i] = np.array([r_norm[j,i], g_norm[j,i], b_norm[j,i]])
    return rgb_grid


##############################################################################
##############################################################################

#%%

######################################PLOTTING RESULTS########################
##############################################################################

def plot_Eu_dict_and_composite(rgb_grid, Eu_sol_dict):
    
    plt.figure()
    ##allowed wavelengths = 410,443,486,551,638,671
    Eu_plot_subplot_loc = [5,6,7,9,10,11] ##subplot location list 
    i = 0 ##index for the subplot location 
    for lam in wavelengths:
        plt.subplot(3,4,Eu_plot_subplot_loc[i]) ##uses the location as given by above
        i += 1
        color_map = plt.cm.get_cmap('hsv')
        color_map.set_under('black')
        if lam == 443 :
            color_map.set_over('white')
            plt.pcolormesh(Eu_sol_dict[lam], cmap=color_map,vmin=0,vmax=.06)  
        else: 
            plt.pcolormesh(Eu_sol_dict[lam], cmap=color_map,vmin=0) 
        plt.colorbar()
        plt.title(r'Eu at Surface ($\lambda$ = {})'.format(lam))
        ##plotting the absorbtion and scattering at each wavelength text 
        ab_wat = abscat(lam,'water')
        ab_diat = abscat(lam, 'Diat')
        ab_syn = abscat(lam, 'Syn')
        x = 110
        y = 150
        text_separation_y = 12
        font_size = 10
        plt.text(x,y,r'$a_{wat} = $ %.4f'%ab_wat[0], fontsize = font_size)
        y = y - text_separation_y
        plt.text(x,y, r'$b_{wat} = $ %.4f'%ab_wat[1], fontsize = font_size)
        y = y - text_separation_y
        plt.text(x,y, r'$a_{diat} = $ %.4f'%ab_diat[0], fontsize = font_size)
        y = y - text_separation_y
        plt.text(x,y, r'$b_{diat} = $ %.4f'%ab_diat[1], fontsize = font_size)
        y = y - text_separation_y
        plt.text(x,y, r'$a_{syn} = $ %.4f'%ab_syn[0], fontsize = font_size)
        y = y - text_separation_y
        plt.text(x,y, r'$b_{syn} = $ %.4f'%ab_syn[1], fontsize = font_size)
    
        # if i > 5: ##The axis will only be on the first graph 
        plt.axis('off')
    
    ##The rainbow with markers of wavelengths used 
    plt.subplot(3,10,(10,30)) ##making it a small top thing at the top of the array
    NP = 750 - 380
    P = np.linspace(381,750,NP)
    P_rgb = np.zeros((NP,3))
    for k in range(NP): 
        P_rgb[k] = wavelength_to_rgb.wavelength_to_rgb(P[k])
    
    len_x = 100
    len_y = len(P)
    T = np.zeros((len_y,len_x,3)) ##slightly confusing but matrix is oriented this way 
    for i in range(len_x):
        T[:,i,:] = P_rgb / 255
    
    plt.imshow(T,origin = 'lower', aspect=('auto'))
    # plt.axis('off')
    indexes = []
    for lam in wavelengths: 
        index = lam - 380
        plt.plot([0,len_x-1], [index,index], 'k')
        # plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
        indexes.append(index)
    plt.yticks(indexes, wavelengths)
    plt.xticks([])
    
    ##Plotting the composite rgb_grid image 
    plt.subplot(3,5,2)
    plt.imshow(rgb_grid, origin=('lower'), aspect=('auto'))
    plt.title('True Color \n Wavelength Composite')
    plt.axis('off')
    # plt.colorbar()
    
    ##plotting the diat at surface concen. and syn at surf. conc. 
    ##gettin gth eunits 
    diat_units = roms_data.variables['diatom'].units
    nano_units = roms_data.variables['nanophytoplankton'].units
    diatom = roms_data.variables['diatom'][0,-1,:,:] 
    nanophyt = roms_data.variables['nanophytoplankton'][0,-1,:,:] 
    ##from Nitogen to Chl
    chl_diatom = (diatom*Chl2NL)  ## Miles Will confirm this.  
    chl_nanophyt = (nanophyt*Chl2NL)
    
    plt.subplot(3,5,1)
    OCx_sol = ocean_color_sol(Eu_sol_dict)
    color_map = plt.cm.get_cmap('hsv')
    plt.pcolormesh(OCx_sol, cmap = color_map.reversed(), norm=mpl.colors.LogNorm(), vmin=.01, vmax = 20)
    plt.title('Nasa Ocean Color(OCx_sol)')
    plt.colorbar()
    plt.axis('off' )
    
    plt.subplot(3,5,3)
    plt.pcolormesh(diatom, norm = mpl.colors.LogNorm())
    plt.title('Log Plot Diatom Surface Concentration \n [{}]'.format(diat_units))
    plt.colorbar()
    plt.axis('off')
    
    plt.subplot(3,5,4)
    plt.pcolormesh(nanophyt, norm = mpl.colors.LogNorm())
    plt.title('Log Plot Nanophytoplankton(Syn) \n Surface Concentration \n [{}]'.format(nano_units))
    plt.colorbar()
    plt.axis('off')
    
    ##plotting some date text. The date from ocean_time. 
    plt.figtext(.02,.95,'ROMS Date and Time: \n {}'.format(roms_date[0]))
    ##plotting some units text for absorbtion and scattering
    plt.figtext(.02,.04, 'a_wat [m^-1] \nb_wat [m^-1] \na_phy [m^2 / mg Chl] \nb_phy [m^2 / mg Chl]')
    
    return

# plot_Eu_dict_and_composite(rgb_grid, Eu_sol_dict)
##############################################################################
##############################################################################

#%%
#############Nasa Ocean Color Study###########################################
##############################################################################
##wavelengths = [410,443,486,551,638,671] ##Taken from John Wilkin email regarding sattelite waveelengths
def ocean_color_sol(Eu_sol_dict): 
    ##the closest our wavelenghts are to the values of ocean color lam
    ## 443(ours) = 443(theirs); 551(ours) = 555(theirs); 671(ours) = 670(theirs) 
    lam_b = 443
    lam_g = 551
    lam_r = 671
    
    ##getting the R_rs from upwelling and the initial conditiongs
    
    R_rs_b = ocean_color.R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_b]) ##R_rs as a function of the blue wavelength 
    R_rs_g = ocean_color.R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_g]) ##green 
    R_rs_r = ocean_color.R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_r]) ##red 
    
    
    OCx_sol = ocean_color.OCx_alg(R_rs_b, R_rs_g, lam_b, lam_g)
    CI_alg_sol = ocean_color.CI_alg(R_rs_r, R_rs_g, R_rs_b, lam_r, lam_g, lam_b)
    
    return OCx_sol


def ocean_color_plot(Eu_sol_dict) :
    ##the closest our wavelenghts are to the values of ocean color lam
    ## 443(ours) = 443(theirs); 551(ours) = 555(theirs); 671(ours) = 670(theirs) 
    lam_b = 443
    lam_g = 551
    lam_r = 671
    
    ##getting the R_rs from upwelling and the initial conditiongs
    
    R_rs_b = ocean_color.R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_b]) ##R_rs as a function of the blue wavelength 
    R_rs_g = ocean_color.R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_g]) ##green 
    R_rs_r = ocean_color.R_RS(E_d_0, E_s_0, Eu_sol_dict[lam_r]) ##red 
    
    
    OCx_sol = ocean_color.OCx_alg(R_rs_b, R_rs_g, lam_b, lam_g)
    CI_alg_sol = ocean_color.CI_alg(R_rs_r, R_rs_g, R_rs_b, lam_r, lam_g, lam_b)
    
    plt.figure()
    plt.subplot(121)
    plt.pcolormesh(CI_alg_sol, cmap = 'hsv')
    plt.title('CI_alg_sol')
    plt.colorbar()
    plt.subplot(122)
    color_map = plt.cm.get_cmap('hsv')
    plt.pcolormesh(OCx_sol, cmap = color_map.reversed(), norm=mpl.colors.LogNorm(), vmin=.01, vmax = 20)
    plt.title('OCx_sol')
    plt.colorbar()

    return 

# ocean_color_plot(Eu_sol_dict)

################################################################################
################################################################################

#%%
# total_water_points = 0
# for mask in maskr.flatten() : 
#     if mask == 1: 
#         total_water_points += 1 

# bad_points_gt_perc_dict = {}
# bad_points_lt_perc_dict = {}
# for lam in wavelengths : 
#     bad_points_gt_perc_dict[lam] = bad_points_gt_dict[lam] / total_water_points
#     bad_points_lt_perc_dict[lam] = bad_points_lt_dict[lam] / total_water_points

    

        
#%%

########################################PLOTTING Single ROMS Profile##########
##############################################################################
# def smoothing(val):

#     val_3_avg = np.zeros(len(val))
#     for i in range(1,len(val) - 1): 
#         val_3_avg[i] = (val[i-1] + val[i] + val[i+1]) / 3.0
        
#     ##the boundaries 
#     val_3_avg[0] = val[0]
#     val_3_avg[-1] = val[-1]
#     return val_3_avg


# lam = 443

# ab_wat = abscat(lam,'water')
# ab_diat = abscat(lam, 'Diat')
# ab_syn = abscat(lam, 'Syn')

# # j,i = 177, 90 ##bok point
# # j,i = 179, 17 ##good 
# j,i = 77,90
# # j, i = 160,90
# # j,i = 150,1

# diatom = roms_data.variables['diatom'][0,:,j,i] 
# nanophyt = roms_data.variables['nanophytoplankton'][0,:,j,i] 

# z_r0 = z_r[:,j,i] ##the verticle array of z_r, the cell centers 
# z_w0 =  z_w[:,j,i] ##z_w is the cell edges 

        
# ##from Nitogen to Chl
# chl_diatom = (diatom*Chl2NL) #+ .5 ## Miles Will confirm this.  
# chl_nanophyt = (nanophyt*Chl2NL) #+ .5

# # plt.figure()
# # plt.plot(diatom, z_r0)



# zarr = np.linspace(-200,0,150) ##1m res to 200m 

# if z_w[0,j,i] > zarr[0]:
#     z_w0 = np.linspace(z_w[0,j,i],z_w[-1,j,i], 1*int(z_w[-1,j,i] - z_w[0,j,i]))
#     # z_w0 = np.linspace(z_w[0,j,i],z_w[-1,j,i], int(len(zarr)))
# elif z_w[0,j,i] < zarr[0]:
#     z_w0 = zarr

# ##smoooth
# # chl_diatom = smoothing(chl_diatom)
# # chl_nanophyt = smoothing(chl_nanophyt)




# chl_w_diatom = np.interp(z_w0,z_r0,chl_diatom) ##interpolating the chl to the z_w grid
# chl_w_nanophyt = np.interp(z_w0,z_r0,chl_nanophyt)
# phy_prof_and_ab_tuple_list = [(chl_w_diatom,ab_diat[0], ab_diat[1]), 
#                               (chl_w_nanophyt, ab_syn[0], ab_syn[1])]
# # ROMS_ocean_irr = irr.ocean_irradiance(z_w0,z_w0, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)
# # [(0*z_w0,0,0)]
# ROMS_ocean_irr = irr.ocean_irradiance(z_w0,z_w0, ab_wat,phy_prof_and_ab_tuple_list, E_d_0, E_s_0)
# Ed, Es, Eu = ROMS_ocean_irr

# plt.figure()

# plt.subplot(131)
# plt.plot(chl_w_diatom, z_w0,label='diatom')
# plt.plot(chl_w_nanophyt, z_w0, label='nano')
# plt.title("Phyt Profile(wavelength = {})".format(lam))
# plt.grid()
# plt.legend()

# plt.subplot(132)
# plt.title('num_opt')
# plt.plot(Ed, z_w0[:], label='Ed')
# plt.plot(Es, z_w0[1:], label='Es')
# plt.plot(Eu, z_w0[1:], label='Eu')
# plt.grid()
# plt.legend()

# ROMS_ocean_irr_bvp = irr_bvp.ocean_irradiance(z_w0, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)

# Ed, Es, Eu = ROMS_ocean_irr_bvp

# plt.subplot(133)
# plt.title('bvp ')
# plt.plot(Ed, z_w0[:], label='Ed')
# plt.plot(Es, z_w0[:], label='Es')
# plt.plot(Eu, z_w0[:], label='Eu')
# plt.grid()
# plt.legend()

# ##############################################################################
# ##############################################################################

#%%

##########################MAIN################################################
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

if __name__ == "__main__" : 
    Eu_sol_dict, bad_points_gt_dict, bad_points_lt_dict,good_points = make_Eu_dict()
    # Eu_sol_dict = make_Eu_dict()
    rgb_grid = wavelength_composite(Eu_sol_dict)
    plot_Eu_dict_and_composite(rgb_grid, Eu_sol_dict)
    # ocean_color_plot(Eu_sol_dict)
    
    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##############################################################################



#################################################################################
# ##saving some data to CSV for chris 
# a = ab_wat[0] + ab_diat[0] * diatom + ab_syn[0] * nanophyt
# b = ab_wat[1] + ab_diat[1] * diatom + ab_syn[1] * nanophyt

# df = pd.DataFrame({'z_roms' : z_r[:,j,i], "diatom" : diatom, "nanophytoplankton" : nanophyt, "a" : a, 'b' : b})
# df.to_csv("ChrisData.csv", index=False)
################################################################################


###############################################################################
#############################CONDITION NUMBER STUDY############################
# lam = 410

# ab_wat = abscat(lam,'water')
# ab_diat = abscat(lam, 'Diat')
# ab_syn = abscat(lam, 'Syn')

# diatom = roms_data.variables['diatom'][0,:,:,:] 
# nanophyt = roms_data.variables['nanophytoplankton'][0,:,:,:] 

# ##from nitrogen to chl
# chl_diatom = (diatom*Chl2NL)  ## Miles Will confirm this.  
# chl_nanophyt = (nanophyt*Chl2NL)

# cond_numbs = np.zeros(np.size(maskr))
# for j in range(nyi): 
#     for i in range(nxi): 
#          chl_diatom0 = diatom[:,j,i]
#          chl_nanophyt0 = nanophyt[:,j,i]
          
#          z_r0 = z_r[:,j,i]  ##rho point s
#          if z_w[0,j,i] > zarr[0]:
#              z_w0 = np.linspace(z_w[0,j,i],z_w[-1,j,i], 1*int(z_w[-1,j,i] - z_w[0,j,i]))
#          elif z_w[0,j,i] < zarr[0]:
#              z_w0 = zarr
#            # z_w0 = z_w[:,j,i]
#          chl_w_diatom0 = np.interp(z_w0,z_r0,chl_diatom0) ##interpolating the chl to the z_w grid
#          chl_w_nanophyt0 = np.interp(z_w0,z_r0,chl_nanophyt0)
          
          
#          phy_prof_and_ab_tuple_list = [(chl_w_diatom0,ab_diat[0], ab_diat[1]), 
#                                               (chl_w_nanophyt0, ab_syn[0], ab_syn[1])]
          
#          ROMS_ocean_irr_and_cond = irr.ocean_irradiance(z_w0, z_w0, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)
         
#          cond_numbs[i*j] = ROMS_ocean_irr_and_cond[-1]
         
# #%%

# for i in range(len(cond_numbs)) :
#     if cond_numbs[i] == 0 :
#         cond_numbs = np.delete(cond_numbs, i)

# plt.figure()
# plt.hist(cond_numbs)
###############################EXTRA STUFF######3#############################

         ##constructing a z array that is in between the cells of z_w with a constant number cell point in each cell 
# z_cell_res = 100 ##I want four points between each z_w0 edge, which means 5 spaces 
# z = np.zeros((len(z_w0)-1) * (z_cell_res-1)) ##should be the length so that each cell is split into 4 points
# zarr_dz_z = np.zeros(len(z_w0) - 1) ##the dz between the two edges divided by the res
# for k in range(len(zarr_dz_z)): ##zarr_dz_z is different for each cell becasue of roms stretcing 
#       zarr_dz_z[k] = (z_w0[k+1] - z_w0[k]) / z_cell_res 

# i = 0 
# for j in range(len(z_w0[:-1])) : ##starting at zbot and going to the edge before ztop. 
#     for k in range(1,z_cell_res) : 
#         z[i] = z_w0[j] + k*zarr_dz_z[j] 
#         i +=1 
 

##STRETCHING STUFF 
# plt.figure()
# plt.plot(diatom, z_r[:,j,i],label='diatom')
# plt.plot(nanophyt, z_r[:,j,i], label='nano')
# plt.legend()

# def stretching(z, dz_min, dz_max,k) : 
#     """
#     z is a negative scalar point from which dz is returned
#     """
#     alpha = dz_max - dz_min 
    
#     dz = dz_min + (alpha * -z) / (k + -z)
    
#     return dz 

# k = 100 ##m half value
# dz_min = .001
# dz_max = 10

# zarr_stretched = []

# ztop = 0 #z_w0[-1]
# zbot = z_w0[0]

# zarr_stretched.append(ztop)
# while True :
#     z = zarr_stretched[-1]
#     new_z = z - stretching(z, dz_min, dz_max,k)
#     print(z)
#     if new_z > zbot : 
#         zarr_stretched.append(new_z)
#     elif new_z < zbot : 
#         break 



         

            

