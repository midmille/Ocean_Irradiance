# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 07:44:02 2021

@author: miles

ABOUT
----- 
This file is to realize the sensitivity of R_rs and chl_a as calculated from 
OCx on the initial surface values of Ed and Es... So Ed0 and Es0. 


"""

## User modules [set PYTHONPATH var]
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
import ocean_irradiance_module.Read_ROMS_Out as RRO
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat

## External modules
import numpy as np
import matplotlib.pyplot as plt 


def Plot_R_rs_vs_Ed0_Es0_Study(Ed0s, Es0s, lam_b, lam_g, R_rs_b, R_rs_g, R_rs_ratio, 
                               chl_OCx, diatom_conc, nano_conc):
    
    fig, axes = plt.subplots(2,1, sharex=True)
    
    ## R_rs Plot
    ## -----------
    im1 = axes[0].plot(Ed0s, R_rs_b, 'b', label=f'R_rs for Blue {lam_b} [nm]')
    im2 = axes[0].plot(Ed0s, R_rs_g, 'g', label=f'R_rs for Green {lam_g} [nm]')

    ## This is the twin x axes, It will have a different y coordinate, mainly the ratio of R_Rs.
    axt_0 = axes[0].twinx()
    im3 = axt_0.plot(Es0s, R_rs_ratio, 'k', label=r'$\frac{\mathrm{R_{rs} \: Blue}}{\mathrm{R_{rs} \: Green}}$')
    
    ## This is the twin y axes. 
    ## This is just so that Es0s will be the upper x-coordinate on the top plot.
    axt_1 = axes[0].twiny()
    ## This come from the assterion that Es0 and Ed0 must sum to 1. 
    axt_1.set_xticks(1.0 - axes[0].get_xticks())
    ## Tightening the layout and correcting for the tick direction.
    axt_1.set_xlim(Es0s.max(), Es0s.min())
    
    
    ## Labels: 
    axes[0].set_title('R_rs Sensitivity')
    axes[0].set_ylabel(r'R_rs [s$r^{-1}$]')
    axt_0.set_ylabel(r'$\frac{\mathrm{R_{rs} \: Blue}}{\mathrm{R_{rs} \: Green}}$')
    axt_1.set_xlabel('Es0', loc='left')
    
    ## Making the legend have all labels
    ims = im1+im2+im3
    im_lbs = [im.get_label() for im in ims]
    axes[0].legend(ims, im_lbs, loc=0)
    
    ## Grid
    axes[0].grid()

    
    
    ## Chl_OCx Plot 
    ## -------------
    axes[1].plot(Ed0s, chl_OCx,'k')
    axes[1].set_title('Chl_OCx')
    axes[1].set_xlabel('Ed0', loc='left')
    axes[1].set_ylabel(r'Chl_a [mg $m^-3$]')
    axes[1].set_xticks(axes[0].get_xticks())
    axes[1].set_xlim(Ed0s.min(), Ed0s.max())
    axes[1].grid()
    
    ## Location of the text at about lower 20% of axes instance.
    x_txt = Ed0s.min() + .1 * (Ed0s.max() - Ed0s.min())
    y_txt = chl_OCx.min() + .1 * (chl_OCx.max() - chl_OCx.min())
    axes[1].text(x_txt, y_txt, f' Constant Diatom [{diatom_conc} mg chl m^-3] \n Constant Nanophyt [{nano_conc} mg chl m^-3]')
    
    return 

    
def R_rs_vs_Ed0_Es0_Study(N, diatom_prof, nano_prof, wavelengths, N_ind, z_phy, 
                          Ed0s, Es0s, Euh, lam_b, lam_g, plot=False): 
    """
    This function performs a study displaying the relationship between R_rs and 
    Es0, Ed0. It assumes a constant phytoplankton profile. 

    Parameters
    ----------
    N : Integer
        The number of levels in the phytoplankton concentration profile. 
    diatom_conc : Float
        The constant diatom concentration.
    nano_conc : Float
        The constant nanophytoplankton concentration.
    wavelengths : List
        COntains the required blue and green wavelengths for OCx algorithim.
    N_ind : Integer
        The number of independent varaibles... Ed0s, Es0s. 
    z_phy : 1-D Array 
        The vertical grid of phytoplankton. Does not matter much since constant 
        phytoplankton. 
    Ed0s : 1-D Array [N_ind]
        The independent variable of Ed0. Downward direct at surface. 
    Es0s : 1-D Array [N_ind]
        The independent variable of Es0. Downward diffuse at surface. Algorithim assumes
        that Es0s = 1 - Ed0s.
    Euh : Float
        Upwelling irradiance bottom boundary condition. Normally taken to be zero.
    lam_b : Float
        The blue wavelength.
    lam_g : Float
        The green wavelength. 
    plot : Boolean, optional
        Flag for plotting result or not. The default is False.

    Returns
    -------
    R_rs_b : 1-D Array [N_ind]
        Remote sensed reflectance for blue wavelength.
    R_rs_g : 1-D Array [N_ind]
        Remote sensed reflectance for green wavelength. 
    R_rs_ratio : 1-D Array [N_ind]
        Remote sensed reflectance ratio. [R_rs Blue] / [R_rs Green]
    chl_OCx : 1-D Array [N_ind]
        Surface chloropyll concentration as calculated from R_rs using the OCx algorithim. 

    """
    
    ## Two species of phytoplankton to be included in Phy object.
    phy_conc = np.zeros((N,2))
    phy_conc[:,0] = diatom_prof
    phy_conc[:,1] = nano_prof

    ## The empty dictionary such that wavelengths will be keys.
    Eu_surf_dict = {}
    
    for lam in wavelengths:
        ## Initiazation of Eu array of equal length as Ed0s, Es0s.
        Eu_surf_dict[lam] = np.zeros(N_ind)
        
        ## The wavelength dependent absorbotion and scattering coefficients
        ab_wat = abscat(lam, 'water')
        ab_diat = abscat(lam, 'Diat') 
        ab_syn = abscat(lam, 'Syn')
        
        ## Creating the Phy object.
        phy = OI.Phy(z_phy, phy_conc, [ab_diat[0], ab_syn[0]], [ab_diat[1], ab_syn[1]] )
        
        for k, (Ed0, Es0) in enumerate(zip(Ed0s, Es0s)) :
            ## Calculating Eu
            Eu_surf_dict[lam][k] = OI.ocean_irradiance(z_phy[0], Ed0, Es0, Euh, ab_wat, 
                                                    phy = phy, N = 30, 
                                                    pt1_perc_zbot = True)[2][-1]
    
    ## Calcualating the respective R_rs values.
    R_rs_b = np.zeros(N_ind)
    R_rs_g = np.zeros(N_ind)
    for k, (Ed0, Es0) in enumerate(zip(Ed0s, Es0s)):
        R_rs_b[k] = OIR.R_RS(Ed0, Es0, Eu_surf_dict[lam_b][k])
        R_rs_g[k] = OIR.R_RS(Ed0, Es0, Eu_surf_dict[lam_g][k])
        
    ## Ratio of R_rs_b / R_rs_g 
    R_rs_ratio = R_rs_b / R_rs_g
    
    ## OCx chlorophyll concentrations. 
    chl_OCx = OIR.OCx_alg(R_rs_b, R_rs_g, lam_b, lam_g)
    
    ## Plotting 
    if plot: 
        Plot_R_rs_vs_Ed0_Es0_Study(Ed0s, Es0s, lam_b, lam_g, R_rs_b, R_rs_g, 
                                   R_rs_ratio, chl_OCx, diatom_prof[0], nano_prof[0])
    
    return R_rs_b, R_rs_g, R_rs_ratio, chl_OCx


def Plot_Chl_OCx_V_Phy_Conc_Ed0_Study(Ed0s, phy_ind, chl_OCx, rel_chl_OCx):
    


    fig, axes = plt.subplots(1,2, sharey=(True))

    ## Plot the true Chl_OCx.
    im1 = axes[0].pcolormesh(Ed0s, phy_ind, chl_OCx, shading = 'auto')
    fig.colorbar(im1, ax = axes[0], label=r'Chl_OCx [mg chl_a $m^{-3}$]')
    axes[0].set_title('Chl_OCx as a Function of \n Phytoplankton Concentration,' +
                      ' Ed0, and Es0')
    axes[0].set_xlabel('Ed0')
    axes[0].set_ylabel(r'Phy Concentration [mg chl_a m^{-3}]')
    
        
    ## Plot relative chl_OCx to Ed0=1 later.
    im2 = axes[1].pcolormesh(Ed0s, phy_ind, rel_chl_OCx, shading = 'auto')
    fig.colorbar(im2, ax = axes[1], label=r'$\frac{\mathrm{chl_OCx - chl_OCx(Ed0=max)}}{\mathrm{chl_OCx(Ed0=max)}}$')
    axes[1].set_title('Relative Chl_OCx as a Function of \n Phytoplankton Concentration,' +
                      ' Ed0, and Es0')
    # axes[1].set_xlabel('Ed0')
    axes[1].set_ylabel(r'Phy Concentration [mg chl_a $m^{-3}$]')
    
    return 


def Chl_OCx_V_Phy_Conc_Ed0_Study(diatom_concs, nano_concs, N_ind, N, wavelengths,
                           z_phy, Ed0s, Es0s, Euh, lam_b, lam_g, diatom_or_nano = 'diatom',
                           plot=True):
      
    if diatom_or_nano == 'diatom':
        phy_ind = diatom_concs
    elif diatom_or_nano == 'nano':
        phy_ind = nano_concs
    
    ## 2D Array of chl_OCx to be stored in.
    ## Phytoplankton concentration will be the vertical row coordinate.
    ## Ed0, Es0 will be the horizontal column coordinate.
    chl_OCx = np.zeros((N_ind,N_ind))
    
    for k, (diatom_conc, nano_conc) in enumerate(zip(diatom_concs, nano_concs)): 
        
        ## Creating the constant profiles. 
        diatom_prof = np.full(N, diatom_conc)
        nano_prof = np.full(N, nano_conc)
        
        ## Calculating Chl_OCx
        chl_OCx[k,:] = R_rs_vs_Ed0_Es0_Study(N, diatom_prof, nano_prof, wavelengths, 
                                             N_ind, z_phy, Ed0s, Es0s, Euh, lam_b, 
                                             lam_g, plot=False)[3]
    
    ## Calculate the relative Chl_OCx from Chl_OCx of max Ed0. 
    rel_chl_OCx = np.zeros_like(chl_OCx)
    for k, Ed0 in enumerate(Ed0s):
        ## Rows represent each concentration level
        ## The 'truth ' is taken to be the Chl_OCx value at the maximum Ed0 value
        ## This could be 100% Ed0 but is not right now because of the divide by zero error.
        ## !! NOTE The max Ed0 is assumed to be the last in the array !!
        rel_chl_OCx[:,k] = (chl_OCx[:,k] - chl_OCx[:,-1] ) / chl_OCx[:,-1]
    
        
    ## Plotting result.
    if plot: 
        Plot_Chl_OCx_V_Phy_Conc_Ed0_Study(Ed0s, phy_ind, chl_OCx, rel_chl_OCx)
    
    
    return     


if __name__ == '__main__':
    
    
    # import argparse
    # parser = argparse.ArgumentParser(description='R_rs sensitivity study.')
    # parser.add_argument('file', help = "Complete Path to ROMS nc file" )
    # parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    # args = parser.parse_args()
    
    ## Experiment Flags
    R_rs_flag = True
    diatom_flag = True
    nanophyt_flag = True
    

    ## The arrays of Ed0 and Es0, mutually dependent sum to 1.
    N_ind = 40
    Ed0s = np.linspace(.1,.9,N_ind)
    Es0s = 1 - Ed0s 
    Euh = 0
    
    ## Blue and Green wavelengths.
    lam_b = 443
    lam_g = 551
    wavelengths = [lam_b,lam_g]
    

    ## The number of levels in profile [same as west coast domain]
    N = 42
    ## The vertical grid that the concentrations are taken.
    z_phy = np.linspace(-1000, 0, N)
    
    
    if R_rs_flag: 
        ## The constant diatom/nanophyt profiles 
        diatom_conc = .5
        diatom_prof = np.full(N, diatom_conc)
        nano_conc = .25 
        nano_prof = np.full(N, nano_conc)
    
        
        R_rs_vs_Ed0_Es0_Study(N, diatom_prof, nano_prof, wavelengths, N_ind, z_phy, 
                              Ed0s, Es0s, Euh, lam_b, lam_g, plot=True)
    
    ## Testing chl_OCx as a function of both phytoplankton concentration and Ed0,Es0. 

    if diatom_flag: 
        
        diatom_or_nano = 'diatom'
        
        ## Second Independent variables is diatom concentration
        diatom_concs = np.linspace(.1,10,N_ind)
        ## Nanophyt is ignored.
        nano_concs = 0 * diatom_concs
        
        ## Running the study
        Chl_OCx_V_Phy_Conc_Ed0_Study(diatom_concs, nano_concs, N_ind, N, wavelengths,
                           z_phy, Ed0s, Es0s, Euh, lam_b, lam_g, diatom_or_nano = diatom_or_nano ,
                           plot=True)
        
    if nanophyt_flag: 
        
        diatom_or_nano = 'nano'
        
        ## Second independent variable is nanophytoplankton concentration
        nano_concs = np.linspace(.1,1,N_ind)
        ## Diatom ignored in this case
        diatom_concs = 0 * nano_concs
        
        ## Running the study
        Chl_OCx_V_Phy_Conc_Ed0_Study(diatom_concs, nano_concs, N_ind, N, wavelengths,
                           z_phy, Ed0s, Es0s, Euh, lam_b, lam_g, diatom_or_nano = diatom_or_nano,
                           plot=True)
        
    

    

    
    
    
        
    
    
    


    
    
