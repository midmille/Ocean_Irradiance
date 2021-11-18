# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 05:40:39 2020

@author: Miles Miller
"""

"""
ROMS Grid shooting method solution
"""
from netCDF4 import Dataset 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
import numpy as np
import wavelength_to_rgb
import ROMS_grid
import os
os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share"
import seapy 
import time
import multiprocessing
import nasa_ocean_color_algorithim as ocean_color
import ocean_irradiance_bvp as irr_bvp
import ocean_irradiance_shoot_method as irr_shoot
import datetime as dt 
import time 
import Metric_Study
import random 
import pickle

#####################################INITIAL PARAMS###########################
##############################################################################
file = '/Users/Miles Miller/da_fwd_002.nc'
##BCs
E_d_0 = .7
E_s_0 = 1 - E_d_0
E_u_h = 0 

##possible wavelengths
wavelengths = [410,443,486,551,638,671] ##Taken from John Wilkin email regarding sattelite waveelengths
# wavelengths = [443]
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


##############################################################################
##############################################################################

###################################Global Funcitons###########################
##############################################################################

def analytical_Ed(zarr,c, E_d_0): 
    """
    Parameters
    ----------
    zarr : vertical 1-D array
        The vertical array of the depth from 0 --> negative depth.

    Returns
    -------
    Downward Direct Irradiance
    

    """
   
    Ed = E_d_0*np.exp(c*zarr)
    return Ed 

def Zbot_func(lam, E_d_0):
    zbots = np.linspace(0, -1000, 2001) 
    a_wat, b_wat = abscat(lam, 'water')
    v_d = .9
    c = (a_wat + b_wat) / v_d
    Ed = analytical_Ed(zbots, c, E_d_0)

    for k, Ed in enumerate(Ed) :
        EdoE_d_0 = Ed / E_d_0
        if EdoE_d_0 < .001 :
            stable_zbot = zbots[k]
            return stable_zbot
        
def Make_Stable_zbot_Dict() :

    """
    Makes a stable zbot dictionary for wavelength list 
    --Water only column 
    --finds the point in column such that Ed = .001*E_d_0
    

    Returns
    -------
    stable_zbot_dict [DICT]
        A dictionary for the list of wavelengths. Each wavelength is a key and the scalar
        stable zbot is its value.

    """
            
##1000 m zbot with 1m res, starts at small zbot and goes to large 
    stable_zbot_dict = {}

    for lam in wavelengths :
        stable_zbot_dict[lam] = Zbot_func(lam, E_d_0)
    
    return stable_zbot_dict


def Eu_at_surface_shooting_method(wavelength,Nlayers = 30):
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
    ##couple init params that really should be defned outside func
    ##however for the time being they are defined within because of the 
    ##necessity of the "wavelength" iterable piece of the function for || computing.
    
    def Log_Trans(zbot,Nlayers):
        zarr_exp = np.linspace(np.log(abs(zbot)),0, Nlayers)
        zarr = -np.exp(zarr_exp)
        return zarr
    
    stable_zbot_dict = Make_Stable_zbot_Dict()
    
    
    ab_wat = abscat(wavelength,'water') ##abscat is the file of all the different 
    ab_diat = abscat(wavelength, 'Diat') ##eyeballed absorbtion and scattering from Dut. 2015
    ab_syn = abscat(wavelength, 'Syn')
    zbot = stable_zbot_dict[wavelength]
    Eu_arr = np.zeros((nyi,nxi)) * np.nan ##Not a Number array so the land is Nan, not 0, helps make land white in pcolormesh
    # zarr = np.linspace(zbot,0,Nlayers) ##10m res to 500m 
    zarr = Log_Trans(zbot,Nlayers)
    diatom = roms_data.variables['diatom'][0,:,:,:] 
    nanophyt = roms_data.variables['nanophytoplankton'][0,:,:,:]
    for j in range(nyi): 
        for i in range(nxi):
            if maskr[j,i] == 1: ##land is a zero, only computes it for water
                 
                ##from Nitogen to Chl
                chl_diatom = (diatom[:,j,i]*Chl2NL)  ## Miles Will confirm this.  
                chl_nanophyt = (nanophyt[:,j,i]*Chl2NS)
                z_r0 = z_r[:,j,i] ##the verticle array of z_r, z_w
                if z_w[0,j,i] > zarr[0]: ##if z_w is shallower than zarr 
                    z_w0 = np.linspace(z_w[0,j,i],0, Nlayers) ##10m Res
                    # z_w0 = np.linspace(z_w[0,j,i],z_w[-1,j,i], int(len(zarr)))
                elif z_w[0,j,i] < zarr[0]: ##if z_w is deeper than zarr, 
                    z_w0 = zarr
                chl_w_diatom = np.interp(z_w0,z_r0,chl_diatom) ##interpolating the chl to the z_w grid
                chl_w_nanophyt = np.interp(z_w0,z_r0,chl_nanophyt)
                phy_prof_and_ab_tuple_list = [(chl_w_diatom,ab_diat[0], ab_diat[1]), 
                                              (chl_w_nanophyt, ab_syn[0], ab_syn[1])]
                N = len(z_w0)
                Nm1 = N-1
                ##initial guess all 
                init_guess = .3
                Ed1 = np.full(N, init_guess)
                Es1 = np.full(N, init_guess)
                Eu1 = np.full(N, init_guess) 
                ##must have the correct BC values in 
                Ed1[Nm1] = E_d_0 ##surface 
                Es1[Nm1] = E_s_0 ##surface 
                Eu1[0] = E_u_h ##bottom
                
                Eu_arr[j,i] = irr_shoot.ocean_irradiance(z_w0, Ed1, Es1, Eu1, ab_wat,phy_prof_and_ab_tuple_list)[-1][-1]
                  
    return Eu_arr



################################post processinng and visualization############
#---post_processing---#


def plot_Eu_dict(rgb_grid, Eu_sol_dict):
    
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


####MAIN#####
#------------------------------------------------------------------------------
if __name__ == '__main__':

    ##multi core python implementaiton 
    start = time.time()
    pool = multiprocessing.Pool(processes = 3)
    Eu_sol = (pool.map(Eu_at_surface_shooting_method, wavelengths))
    end = time.time()
    time_elapsed = start  - end 
    print('runtime:', time_elapsed)
    Eu_sol_dict_shoot_RK4_30Layer_LogTrans = {}
    for k in range(len(Eu_sol)) :
        Eu_sol_dict_shoot_RK4_30Layer_LogTrans[wavelengths[k]] = Eu_sol[k]
    ##LAST RUN TIME: 21 minutes 
    pickle.dump(Eu_sol_dict_shoot_RK4_30Layer_LogTrans, open("Eu_sol_dict_shoot_RK4_30Layer_LogTrans", "wb"))
    rgb_grid = wavelength_composite(Eu_sol_dict_shoot_RK4_30Layer_LogTrans)

    plot_Eu_dict(rgb_grid, Eu_sol_dict_shoot_RK4_30Layer_LogTrans)

#    ocean_color_plot(Eu_sol_dict_shoot_RK4_200Layer_pt1perc_lightLevel)



























