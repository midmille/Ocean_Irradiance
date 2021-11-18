# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 06:14:01 2021

@author: miles

This file is for the quantification of the dependence of the irradiance solution upon 
the choice in depth that the bvp solution iterates to. 

--> The study will be completed for pure water and for an artificial phytoplankton profile

"""


## User Modules
from ocean_irradiance_module import Ocean_Irradiance as OI
from ocean_irradiance_module import Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.PARAMS import Param_Init 
from ocean_irradiance_module import Wavelength_To_RGB as WTRGB

## External Modules
import numpy as np
import matplotlib.pyplot as plt 


if __name__ == '__main__':
    
    ## Initials
    ## ---------
    
    ## Flags 
    wat_or_phy='phy'
    ## independent variable zbot as a fraction of zbot_pt1_perc_light_level
    N_ind = 20
    ## Number of vertical levels 
    N = 30
    ## by appending 1 to the end of the array, the last value in solution will be taken as truth.
    zbot_frac = np.append(np.linspace(.5, 2, N_ind-1), 1)
    
    ## Initializing parameters
    PI = Param_Init()
    PI.wavelengths = [443]
    
    
    ## Study 
    ## ------
    ## Loop wave length
    Eu_array = np.zeros((N,N_ind))
    Eu_dict = {}
    for k, lam in enumerate(PI.wavelengths):
        ## absorbtion/scattering for water. 
        ab_wat = abscat(lam, 'water')
        a_wat, b_wat = ab_wat
        ## finding the .1% light level
        pt1_perc_zbot = OI.zbot_func(PI.Ed0, a_wat, b_wat, PI.v_d)
        ## The different zbots
        zbots = pt1_perc_zbot * zbot_frac
        
        ## Looop of independent variable
        if wat_or_phy == 'wat': 
            for j, zbot in enumerate(zbots):
                Eu_array[:,j] = OI.ocean_irradiance(zbot, PI.Ed0, PI.Es0, PI.Euh, 
                                                 ab_wat, PI.coefficients, phy=None,
                                                 pt1_perc_zbot=False)[2]
            Eu_dict[lam] = Eu_array
            ## The truth is taken to be the Eu solved for using scipy solver.
            Eu_truth = OI.ocean_irradiance(pt1_perc_zbot, PI.Ed0, PI.Es0, PI.Euh, 
                                                 ab_wat, PI.coefficients, phy=None,
                                                 pt1_perc_zbot=False, use_bvp_solver=True)[2]
               
        elif wat_or_phy == 'phy':
            for j, zbot in enumerate(zbots):
                N_z = 100
                z_phy = np.linspace(0, zbot, N_z)
                ## Full profile of 1
                phy_prof = np.full(N_z, 1)
                ab_phy = abscat(lam, 'Diat')
                phy = OI.Phy(z_phy, phy_prof, ab_phy[0], ab_phy[1])
                
                Eu_array[:,j] = OI.ocean_irradiance(zbot, PI.Ed0, PI.Es0, PI.Euh, 
                                                 ab_wat, PI.coefficients, phy=phy,
                                                 pt1_perc_zbot=False)[2]
            Eu_dict[lam] = Eu_array
            ## The truth is taken to be the Eu solved for using scipy solver.
            Eu_truth = OI.ocean_irradiance(pt1_perc_zbot, PI.Ed0, PI.Es0, PI.Euh, 
                                     ab_wat, PI.coefficients, phy=None,
                                     pt1_perc_zbot=False, use_bvp_solver=True)[2]
            


    ## Calculationg the relative difference. 
    
    fig, ax = plt.subplots()
    ax.set_xlabel("Fraction of .01% pure water light level zbot [m]")
    ax.set_ylabel("Max Relative Difference in Eu profile")
    ax.set_title("Eu Profile Sensitivity to Choice of Zbot")
    ax.grid()
    
    rel_diff = np.zeros(N_ind)
    for k, lam in enumerate(PI.wavelengths):
        for j, zbot in enumerate(zbots):
            Eu_truth[0] = Eu_truth[1]
            rel_diff[j] = max(abs(( Eu_dict[lam][:,j] - Eu_truth )))
        
        ax.plot(zbot_frac[:-1], rel_diff[:-1], color=[c/255 for c in WTRGB.wavelength_to_rgb(lam)],
                label=lam)
        
    ax.legend()

    
    