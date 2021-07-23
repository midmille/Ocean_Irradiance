# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 05:52:04 2021

@author: Miles Miller

Zbot Dependence of Three Stream Irradiance Solution
----------------------------------------------------

--> Solves the three stream BVP for an array of points in ROMS output and for an
array of different zbots. 
--> The set of point is chosen from ROMS output using Integrate_and_Bin_Representative_ROMS_Points
module.
--> Calculates the percent difference from SciPy BVP solver. 
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
import Integrate_and_Bin_Representative_ROMS_Points as Rep_Bin
## Importing other modules. 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import numpy as np


##zbot study 
##Normalized version 
# wavelengths = [443]
# zbots =  np.arange(-15, -800,-5)
## a small call to see what the depths of the profiles from rep_points_index 

def True_Eu(): 
    
    return 

def Zbot_Sensitivity_ROMS(rep_points_index, ab_wat, ab_diat, ab_syn, 
                          chl_diatom, chl_nanophyt, mask, z_r, E_d_0, E_s_0, v_d,
                          E_u_h, N, N_zbots) :
    
    ## Creating an empty array for the z grid of given phy profiles to reside in.
    phy_prof_z_arr = np.zeros((len(rep_points_index[:,0]), len(z_r[:,0,0])))
    ## The fraction of zbot at the .1% light level we want to test. 
    fract_zbot = np.linspace(-.5, -5, N_zbots)
    ## Reconstructing the mask to only include the representative ROMS points. 
    rep_mask = np.zeros_like(mask)
    ## Looping over the representative points index. 
    for (yi,xi) in rep_points_index: 
        ## One denotes a point that will be included in calculation.
        rep_mask[yi,xi] = 1
    
    ## Creating a z_r such that hbot == zbot_fract.
    ## This will be used by the Ocean_Irradiance_ROMS module with the pt1_perc_zbot 
    ## flag set to False so that hbot will be used instead of the .1% light level. 
    z_r_zbot = np.copy(z_r)
    ## An empty array to store the result
    Eu_zbot_fract = np.zeros((len(fract_zbot), mask.shape[0], mask.shape[1]))
    ## Loop over the fractional zbots.
    for k, fract in enumerate(fract_zbot): 
        ## Changing hbot to be pt1_perc_zbot * frac.  
        z_r_zbot[0] = Ocean_Irradiance.zbot_func(Ed0, a_wat, b_wat, v_d) * fract
        ## Getting the surface upwelling irradiance for the representative point in the mask. 
        Eu_zbot_fract[k,:,:] = Ocean_Irradiance_ROMS.Eu_at_surface(mask, ab_wat, ab_diat,
                                                            ab_syn, chl_diatom, chl_nanophyt, 
                                                            z_r, E_d_0, E_s_0, E_u_h, N=30,
                                                            pt1_perc_zbot = False)
    ## Now to get the 'true' value of Eu_zbot
        
    ## Different version that loops just over the representative points and not the mask.
    ## This would entail NOT using the Ocean_Irradiance_ROMS module. 
    for k, fract in enumerate(fract_zbot): 
        for k, (yi,xi) in enumerate(rep_points_index): 
            
            zbot = Ocean_Irradiance.zbot_func(Ed0, a_wat, b_wat, v_d) * fract
            
            phy_profs = np.zeros((len(z_r0),2))
            phy_profs[:,0] = chl_diatom
            phy_profs[:,1] = chl_nanophyt
            
            a = np.array([ab_diat[0], ab_syn[0]])
            b = np.array([ab_diat[1], ab_syn[1]])
            
            phy = Ocean_Irradiance.Phy(z_r0, phy_profs, a, b)    
            
            ocean_irr_sol = Ocean_Irradiance.ocean_irradiance(zbot,E_d_0,E_s_0,E_u_h,
                                                            ab_wat,phy, N=N, 
                                                            pt1_perc_zbot = pt1_perc_zbot)
                    
    
        
        
        
        
    
    
    
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
    Eu_arr =  Ocean_Irradiance_ROMS.Eu_at_surface(rep_mask, ab_wat, ab_diat, ab_syn, 
                                                  chl_diatom, chl_nanophyt, z_r, 
                                                  E_d_0, E_s_0, E_u_h, N=N)
    
    
    
    
    
    
    
    ## Looping over the representative points index. 
    for k,(yi,xi) in enumerate(rep_points_index) : 
        ## The phytoplankton z grid at given index.
        phy_prof_z_arr[:,k] = z_r[:, int(yi), int(xi)]
        ## Diatom and nanophyt concentrations at given index. 
        chl_diatom_i = chl_diatom[:,yi,xi]
        chl_nanophyt_i = chl_nanophyt[:,yi,xi]


    

    # stable_zbot_dict = Make_Stable_zbot_Dict()
    # for lam in wavelengths : 
    #     zbots_dict[lam] = stable_zbot_dict[lam] * (-fract_zbot) ##Normalized zbots
    
    Eu_vs_zbot_sol_dict = {}
    Eu_vs_zbot_sol_bvp_dict = {}
    Eu_prof_shoot_zbot400m = {}
    max_layers = 1500
    Eu_prof_sol_bvp = np.empty((max_layers,len(fract_zbot)))
    Eu_prof_sol_bvp[:,:] = np.NaN
    zbot_zarr = np.empty((max_layers,len(fract_zbot)))
    zbot_zarr[:,:] = np.NaN
    Es_prof_sol_bvp = np.empty((max_layers,len(fract_zbot)))
    Es_prof_sol_bvp[:,:] = np.NaN

    for lam in wavelengths : 
        Eu_sol = np.zeros(len(fract_zbot))
        Eu_sol_bvp = np.zeros(len(fract_zbot))
        for k in range(len(fract_zbot)) :
            solution = irr_phy_prof_shoot(lam, zbots_dict[lam][k], nanophyt, diatom, 
                                          phy_z_grid, E_d_0, E_s_0, Nlayers=200,return_BVP=True)
            
            Eu_sol[k] = solution[2][-1]
            Eu_sol_bvp[k] = solution[5][-1]
            zbot_zarr[:len(solution[5]),k] = solution[6]
            Eu_prof_sol_bvp[:len(solution[2]),k] = solution[2]
            Es_prof_sol_bvp[:len(solution[2]),k] = solution[1]
            Ed_prof_sol_bvp[:len(solution[2]),k] = solution[0]
        Eu_vs_zbot_sol_dict[lam] = Eu_sol
        Eu_vs_zbot_sol_bvp_dict[lam] = Eu_sol_bvp
        
    # for lam in wavelengths : 
    #     Eu_prof_shoot_zbot400m[lam] = irr_water_only_shoot(lam, - 20, E_d_0, E_s_0)[2]
    # zarr = irr_water_only_shoot(lam, - 20, E_d_0, E_s_0)[6]

    fig, axes = plt.subplots(1,3)
    (ax,ax1,ax2) = axes 
    for lam in wavelengths : 
        for k,zbot in enumerate(zbots_dict[lam]) :
            ax.plot(Eu_prof_sol_bvp[:,k], zbot_zarr[:,k], label = fract_zbot[k])
            ax.set_title('Eu')
            ax1.plot(Es_prof_sol_bvp[:,k], zbot_zarr[:,k], label = fract_zbot[k])
            ax1.set_title('Es')
            ax2.plot(Ed_prof_sol_bvp[:,k], zbot_zarr[:,k], label = fract_zbot[k])
            ax2.set_title('Ed')
    for axe in axes : 
        axe.legend()
        axe.grid()
    

    
    ##making a list of the x-values = resulting errors
    Eu_diff_rel = np.zeros((len(fract_zbot),len(wavelengths)))
    for k,lam in enumerate(wavelengths)  :
        true_Eu = true_value(lam, stable_zbot_dict[lam])[2][-1]
        # zbot_for_true = 1000
        # true_Eu = irr_water_only_shoot(lam, zbot_for_true)[5][-1]
        # true_Eu = Eu_vs_zbot_sol_bvp_dict[lam][-1]
        Eu_true_Eu_diff = true_Eu - Eu_vs_zbot_sol_bvp_dict[lam]
        Eu_true_Eu_diff_rel = Eu_true_Eu_diff / true_Eu
        Eu_diff_rel[:,k] = (Eu_true_Eu_diff_rel)
        

    # fig, (ax, ax1) = plt.subplots(1,2)
    fig, ax = plt.subplots()
    ax.set_title("a) Eu at Surface Error vs z-bottom Sensitivity Study Using\n Pure Water Absorbtion and Scattering \
    Shooting Method")
    ax.set_ylabel(" Fraction of zbot such that Ed = .001(Ed0) ")
    ax.set_xlabel('Normalized Relative Error from BVP Solver')
    rgb_from_wavelengths = [wavelength_to_rgb.wavelength_to_rgb(k) for k in wavelengths]
    colors = [(k[0]/255, k[1]/255, k[2]/255) for k in rgb_from_wavelengths] ##making a tuple from 0,1
    for k,lam in enumerate(wavelengths) :
        # ax.plot(Eu_diff_rel[:,k], stable_zbot_dict[lam], label = wavelengths[k], color = colors[k])
        ax.plot(Eu_diff_rel[1:,k], fract_zbot[1:], label = wavelengths[k], color = colors[k]) ##nolonger wavelen dep. y-axis
    # x_minor_ticker = mpl.ticker.MultipleLocator(.05)
    # x_major_ticker = mpl.ticker.MultipleLocator(.1)
    # ax.xaxis.set_minor_locator(x_minor_ticker)
    # ax.xaxis.set_major_locator(x_major_ticker)
    y_major_ticker = mpl.ticker.MultipleLocator(.5)
    ax.yaxis.set_major_locator(y_major_ticker)
    # ax.set_xscale('log10')
    ax.grid()
    ax.legend(title='wavelengths[nm]')
    fig.text(.01,.01,'Miles Miller \n{date}'.format(date = dt.date.today()), fontsize = 8)


def main():
    