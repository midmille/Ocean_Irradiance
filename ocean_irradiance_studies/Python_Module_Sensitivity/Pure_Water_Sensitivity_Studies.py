# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 08:11:53 2020

@author: Miles Miller
"""
"""
Pure seawater only sensitivity studies 
"""

import numpy as np
from absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
import ocean_irradiance_bvp as irr_bvp
import ocean_irradiance_shoot_method as irr_shoot
import Metric_Study
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt 
import time 
import wavelength_to_rgb 

################################PARAMS########################################
E_d_0 = .7
E_s_0 = 1 - E_d_0
E_u_h = 0 

wavelengths = [410,443,486,551,638,671] ##Taken from John Wilkin email regarding sattelite waveelengths

##############################################################################

#############################FUNCTIONS########################################

def irr_water_only_shoot(lam, zbot, E_d_0, E_s_0, N_layers=200, shots=3, zbot_pt1=10) : 
    ##we are changing zarr for this sensitivity study 
    # dz = zbot_pt1 / N_layers 
    zarr = np.linspace(zbot, 0, N_layers) ##1m res 
    # N_layers_zarr = round(abs(zbot / dz))
    # zarr = np.linspace(zbot, 0, N_layers_zarr)
    
    
    ab_wat = abscat(lam,'water')
    # print(lam, ab_wat)      
    N = len(zarr)
    Nm1 = N-1
    ##initial guess all ones 
    Ed1 = np.ones(N)
    Es1 = np.ones(N)
    Eu1 = np.ones(N) 
    ##must have the correct BC values in 
    Ed1[Nm1] = E_d_0 ##surface 
    Es1[Nm1] = E_s_0 ##surface 
    Eu1[0] = E_u_h ##bottom 
    
    ##must make a list of zeros
    phy_prof_and_ab_tuple_list = [(zarr*0,0,0)]
    
    ocean_irradiance_output = irr_shoot.ocean_irradiance(zarr, Ed1, Es1, Eu1, ab_wat, phy_prof_and_ab_tuple_list, shots) 
    ocean_irradiance_bvp_output = irr_bvp.ocean_irradiance(zarr, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)
    Ed, Es, Eu = ocean_irradiance_output
    Ed_bvp, Es_bvp, Eu_bvp =  ocean_irradiance_bvp_output
    return Ed, Es, Eu, Ed_bvp, Es_bvp, Eu_bvp, zarr
    # return Ed, Es, Eu, zarr

def true_value(lam, zbot):
    # zbot = -600
    zarr = np.linspace(zbot, 0, 5000) ##1m res 
    
    ab_wat = abscat(lam,'water')
    # print(lam, ab_wat)
    N = len(zarr)
    Nm1 = N-1
    ##initial guess all ones 
    Ed1 = np.ones(N)
    Es1 = np.ones(N)
    Eu1 = np.ones(N) 
    ##must have the correct BC values in 
    Ed1[Nm1] = E_d_0 ##surface 
    Es1[Nm1] = E_s_0 ##surface 
    Eu1[0] = E_u_h ##bottom 
    
    ##must make a list of zeros
    phy_prof_and_ab_tuple_list = [(zarr*0,0,0)]
    
    ocean_irradiance_output = irr_bvp.ocean_irradiance(zarr, ab_wat, phy_prof_and_ab_tuple_list, E_d_0, E_s_0)
    
    Ed, Es, Eu = ocean_irradiance_output
    
    return Ed, Es, Eu

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


def Stable_zbot(lam,zbots, E_d_0):
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
    
    zbots = np.linspace(0, -1000, 2001) ##1000 m zbot with 1m res, starts at small zbot and goes to large 
    stable_zbot_dict = {}

    for lam in wavelengths :
        stable_zbot_dict[lam] = Stable_zbot(lam, zbots, E_d_0)
    
    return stable_zbot_dict
        

##############################################################################


##################################STUDIES#####################################
##############################################################################
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#

#%% Number of Shots Study 
###############################NUMBER OF SHOTS STUDY##########################
#----------------------------------------------------------------------------#

##number of shots study 
##Must load functions: Stable_zbot and analyitical_Ed 
##getting the stable zbots to use for this study
zbots = np.linspace(-5,-1000, 500)
stable_zbot_dict = {}
for lam in wavelengths: 
    stable_zbot = Stable_zbot(lam, zbots, E_d_0)
    stable_zbot_dict[lam] = stable_zbot 
    
##Now to do different number of shots into the water only shoot method 
Nlayers = 10 ##choosing a constant number of layers for each wavelengtht
shots_list = np.arange(2, 6, 1)
shots_result_dict = {} ##for each wavelength result will be BV at surface for Eu...
for lam in wavelengths : ## making another dict of lam and shots result 
    shots_result_arr = np.zeros(len(shots_list))
    for k,shots in enumerate(shots_list) : 
        solution = irr_water_only_shoot(lam, stable_zbot_dict[lam], E_d_0, E_s_0, Nlayers, shots)
        shots_result_arr[k] = solution[2][-1] ##surface value.. what we are shooting towards
    shots_result_dict[lam] = shots_result_arr
##Results of this study: 
    ##I found that changing the number of shots does not change the solution 
    ##for any more than three shots
    ##even for small Nlayers 
    ##thus 3 shots is perfect, no more needed 

#%% 
###############################METRIC STUDY###################################
#----------------------------------------------------------------------------#
##continuation of number of shots study 
##changing initial guesses so as to better understand how the metric changes 
##and to understan why only three shots are needed for linear differential equations 

stable_zbot_dict = Make_Stable_zbot_Dict()  
lam = 443
zbot = stable_zbot_dict[lam]
##truth 
true_Eu0 =  irr_water_only_shoot(lam, zbot, E_d_0, E_s_0)[2][-1]
metric_array = Metric_Study.metric_study(lam, zbot, E_d_0, E_s_0, E_u_h, true_Eu0, N_layers=200)
zarr = metric_array[0]
Eu0s = metric_array[1]
Fmetric = metric_array[2]
Eu_profs = metric_array[3]


##plotting results 
##metric vs Eu0s 
fig, ax = plt.subplots()
ax.plot(Eu0s, Fmetric)
ax.grid()
ax.set_ylabel('Eu at Bottom (Fmetric)')
ax.set_xlabel('Eu at Top')


##plotting Eu profs vs zarr 
fig, ax = plt.subplots()
for k, Eu0 in enumerate(Eu0s) :
    ax.plot(Eu_profs[:,k], zarr, label = round(Eu0,4))
ax.grid()
ax.set_xlim(-.1, .1)
ax.legend(title='Eu at Surface:')
ax.set_ylabel('z [m]')
ax.set_xlabel('Eu')    
ax.set_title('Eu Profiles for Different Eu0s')


    
#%%
###############################Resolution Sensitivity Study###################
#----------------------------------------------------------------------------#

##Resolution Sensitivity study 
##N layer remeains constant for this 
stable_zbot_dict ##we need this 
Nlayers_list = np.arange(5,300, 5)
dz_study_results_dict = {}
dz_study_time_dict = {}
dz_study_zarr_dict = {}
for lam in wavelengths : 
    dz_result_Eu_arr = np.zeros(len(Nlayers_list))
    dz_result_time = np.zeros(len(Nlayers_list))
    dz_result_zarrs = np.empty((np.max(Nlayers_list), len(Nlayers_list)))
    dz_result_zarrs[:] = np.NaN ##empty NAN array so I can fill to varried zarr levels
    for k,Nlayers in enumerate(Nlayers_list) :
        start = time.perf_counter() ##start timer 
        solution = irr_water_only_shoot(lam, stable_zbot_dict[lam], E_d_0, E_s_0, Nlayers)
        dz_result_Eu_arr[k] = solution[2][-1] ##last index in array is the surface Eu
        dz_result_zarrs[:Nlayers,k] = solution[3] ##zarr of length Nlayers
        stop = time.perf_counter() ##end timer 
        timer = stop - start ##time btw 
        dz_result_time[k] = timer
    dz_study_results_dict[lam] = dz_result_Eu_arr
    dz_study_time_dict[lam] = dz_result_time 
    dz_study_zarr_dict[lam] = dz_result_zarrs

##calculating relative error 
##taking the most Nlayers to be truth 
dz_study_Eu_diff = {}
dz_study_Eu_rel_err = {}
for lam in wavelengths : 
    dz_study_Eu_diff[lam] = dz_study_results_dict[lam][-1] - dz_study_results_dict[lam] ##last one in array minus all others
    print(dz_study_Eu_diff[lam])
    dz_study_Eu_rel_err[lam] = dz_study_Eu_diff[lam] / dz_study_results_dict[lam][-1]
#%% ##above is calculations, below plotting

##plotting the results 
fig,ax = plt.subplots()
ax.set_title('Resolution Sensitivity Study')
rgb_from_wavelengths = [wavelength_to_rgb.wavelength_to_rgb(k) for k in wavelengths]
colors = [(k[0]/255, k[1]/255, k[2]/255) for k in rgb_from_wavelengths] ##making a tuple from 0,1
for k,lam in enumerate(wavelengths) :
    # ax.plot(dz_study_time_dict[lam], dz_study_Eu_rel_err[lam], label = wavelengths[k], color = colors[k])
    # ax.plot(dz_study_time_dict[lam], dz_study_Eu_diff[lam], label = wavelengths[k], color = colors[k])
    ax.plot(Nlayers_list, dz_study_Eu_rel_err[lam], label = wavelengths[k], color = colors[k])
    
    # ax.plot(Nlayers_list, dz_study_Eu_diff[lam], label = wavelengths[k], color = colors[k])
# ax.plot(Nlayers_list, np.full((len(Nlayers_list)),.001), linestyle = '--', color='black', label='1% Relative Error Level')
x_1perc_err_loc = 30
x_pt1perc_err_loc = 90
ax.plot(np.full(50,x_1perc_err_loc), np.linspace(0,.3,50), ls='--',lw = 2, color='k',label='1.0% Relative Error')
ax.plot(np.full(50,x_pt1perc_err_loc), np.linspace(0,.3,50), ls=':',lw = 2, color='k',label='0.1% Relative Error')
ax.set_ylabel('Relative Error to Highest Resolution')
ax.set_xlabel('Computational Cost[Number of Layers]')
ax.set_xticks(np.append(ax.get_xticks(),[x_1perc_err_loc, x_pt1perc_err_loc]))
ax.set_xlim(0,300)
ax.set_ylim(-.01,.32)
# ax.text(130,.11, '1.0% Error Layers = {} \n0.1% Error Layers = {}'.format(x_1perc_err_loc, x_pt1perc_err_loc))
# ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(7))
# ax.set_xticks(np.linspace(0, np.max(dz_study_time_dict[410])), 7)
# ax.set_xlim(0,  np.max(dz_study_time_dict[410]))
ax.legend()

# ax1 = ax.twiny()
# ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(7))
# ax1.set_xlim(Nlayers_list[0], Nlayers_list[-1])
# # ax1.set_xticks(np.linspace(0, np.max(Nlayers_list)), 7)

ax.grid()

#%%
###########################zbot dependence sensitivity########################
#----------------------------------------------------------------------------#

##zbot study 
##Normalized version 
# wavelengths = [443]
# zbots =  np.arange(-15, -800,-5)
zbots_dict = {} ##fraction zbots of the stable zbot dict thus wavelength dep.
fract_zbot = np.arange(-.5,-5.25,-.1) ##to be multiplied by the zbot stable 
# fract_zbot = np.asarray([-.25, -.5, -1])
for lam in wavelengths : 
    zbots_dict[lam] = stable_zbot_dict[lam] * (-fract_zbot) ##Normalized zbots

Eu_vs_zbot_sol_dict = {}
Eu_vs_zbot_sol_bvp_dict = {}
Eu_prof_shoot_zbot400m ={}
max_layers = 1500
Eu_prof_sol_bvp = np.empty((max_layers,len(fract_zbot)))
Eu_prof_sol_bvp[:,:] = np.NaN
zbot_zarr = np.empty((max_layers,len(fract_zbot)))
zbot_zarr[:,:] = np.NaN
Es_prof_sol_bvp = np.empty((max_layers,len(fract_zbot)))
Es_prof_sol_bvp[:,:] = np.NaN
Ed_prof_sol_bvp = np.empty((max_layers,len(fract_zbot)))
Ed_prof_sol_bvp[:,:] = np.NaN
for lam in wavelengths : 
    Eu_sol = np.zeros(len(fract_zbot))
    Eu_sol_bvp = np.zeros(len(fract_zbot))
    for k in range(len(fract_zbot)) :
        solution = irr_water_only_shoot(lam, zbots_dict[lam][k], E_d_0, E_s_0, zbot_pt1=stable_zbot_dict[lam])
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
##truth 
#%%
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

#%%
##Error calculations 

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
    
##object oriented attempt at plotting 
#%%
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
# ax.text(.035,-190, r"Rel. Error = $ \frac{(\mathrm{bvp\ sol.\ with\ .2m\ res.\ at\ zbot=1000m}) - (\mathrm{Eu\ shoot\ method\ 1m\ res.\ for\ z-bot})}{(\mathrm{bvp\ sol.\ with\ .2m\ res.\ at\ zbot=1000m})} $", fontsize = 12)


# ##calculating stable zbots
# zbots2 = np.linspace(0, -1000, 1001) ##1000 m zbot with 1m res, starts at small zbot and goes to large 
# stable_zbot_dict = {}
# zarr_zbot_dict ={}
# for lam in wavelengths :
#     stable_zbot_dict[lam] = Stable_zbot(lam, zbots2, E_d_0)
# Eu_sol_prof_stable_zbot = {}
# for lam in wavelengths : 
#     sol_water = irr_water_only_shoot(lam, stable_zbot_dict[lam], E_d_0, E_s_0)
#     Eu_sol_prof_stable_zbot[lam], zarr_zbot_dict[lam] = sol_water[2], sol_water[6]


# ##plotting 
# # fig, ax = plt.subplots()
# ax1.set_ylabel('z[m]')
# ax1.set_xlabel('Eu')
# ax1.set_title('b) Eu Profile for Each Wavelength at zbot Where Ed = .001*E_d_0')
# rgb_from_wavelengths = [wavelength_to_rgb.wavelength_to_rgb(k) for k in wavelengths]
# colors = [(k[0]/255, k[1]/255, k[2]/255) for k in rgb_from_wavelengths] ##making a tuple from 0,1
# for k,lam in enumerate(wavelengths) :
#     ax1.plot(Eu_sol_prof_stable_zbot[lam], zarr_zbot_dict[lam], label = wavelengths[k], color = colors[k])
# ax1.grid()
# ax1.legend(title='wavelengths[nm]')
# # ax1.set_title("Eu Profile at zbot=100m for Different Wavelenghts")
# # ax1.set_ylabel("z [m]")
# # ax1.set_xlabel('Eu')
# # for k,lam in enumerate(wavelengths[0:6]) :
# #     ax1.plot(Eu_prof_shoot_zbot400m[lam], zarr, label = wavelengths[k], color = colors[k])
# # ax1.grid()
# # ax1.legend(title='wavelengths[nm]')

#%%
########################INITIAL ED0 ES0 STUDY#####################################
#----------------------------------------------------------------------------#

##Calculating the depth for each wavelength such that Ed is .1% of surface value
##for water only 


        
E_d_0s = np.linspace(.001,1,7)
E_s_0s = 1 - E_d_0s
zbots = np.linspace(0, -1000, 1001) ##1000 m zbot with 1m res, starts at small zbot and goes to large 
stable_zbot_dict = {}
##how stable zbot changes with E_d_0 and E_s_0
for lam in wavelengths :
    stable_zbot_arr = np.zeros(len(E_d_0s))
    for k,E_d_0 in enumerate(E_d_0s) :
        stable_zbot_arr[k] = Stable_zbot(lam, zbots, E_d_0)
    stable_zbot_dict[lam] = stable_zbot_arr
##how chnages with E_s_0 chnages Eu profiles 
N_layers = 200
Eu_prof_zbot_and_E_d_0s_dict ={}
zarrs_dict = {}
for lam in wavelengths :
    Eu_profiles = np.zeros((len(E_d_0s), N_layers))
    zarrs = np.zeros_like(Eu_profiles)
    for k,zbot in enumerate(stable_zbot_dict[lam]) :
        
        solution = irr_water_only_shoot(lam, zbot, E_d_0s[k], E_s_0s[k], N_layers)
        Eu_profiles[k] = solution[2]
        zarrs[k] = solution[6]
        
    Eu_prof_zbot_and_E_d_0s_dict[lam] = Eu_profiles##Eu result from shoot is third inreturn tuple
    zarrs_dict[lam] = zarrs   
##plotting Eu profiles for each wavelength for each E_s_0 and E_d_0
fig, axes = plt.subplots(2,3)
axes = axes.flatten()
for k,lam in enumerate(wavelengths) : 
    for j in range(len(zarrs_dict[lam])) :
        axes[k].plot(Eu_prof_zbot_and_E_d_0s_dict[lam][j], zarrs_dict[lam][j], label=round(E_s_0s[j], 4))
    axes[k].grid()
    axes[k].set_title(wavelengths[k])
    if k==0 or k==3 :
        axes[k].set_ylabel('z [m]')
    axes[k].set_xlabel('Eu')
axes[k].legend(title='E_s_0')
fig.text(.3,.92,'Eu profiles at stable zbot for variable Es0 and Ed0', fontsize = 18)
    

##plotting result of stable zbot for diff E_d_0s and E_s_0s 
fig, ax = plt.subplots()
rgb_from_wavelengths = [wavelength_to_rgb.wavelength_to_rgb(k) for k in wavelengths]
colors = [(k[0]/255, k[1]/255, k[2]/255) for k in rgb_from_wavelengths] ##making a tuple from 0,1
for k,lam in enumerate(wavelengths) :
    ax.plot(E_d_0s, stable_zbot_dict[lam], label = wavelengths[k], color = colors[k])
ax.grid()
ax.set_xlim(0,1)
ax.legend(title='wavelengths')
ax.set_ylabel('zbot where Ed/E_d_0 = .001 [m]')
ax.set_xlabel('E_d_0')

ax1 = ax.twiny()
ax1.set_xlabel('E_s_0')
ax1.set_xlim(1,0)

ax.set_title('zbot vs E_d_0 and E_s_0', loc='left', fontsize=14)

        


#%%
##how Eu at surface changes with wavelength dependent zbot
zbots = np.linspace(0, -1000, 1001) ##1000 m zbot with 1m res, starts at small zbot and goes to large 
stable_zbot_dict = {}
zarr_zbot_dict ={}
for lam in wavelengths :
    stable_zbot_dict[lam] = Stable_zbot(lam, zbots, E_d_0)
Eu_sol_prof_stable_zbot = {}
for lam in wavelengths : 
    sol_water = irr_water_only_shoot(lam, stable_zbot_dict[lam], E_d_0, E_s_0)
    Eu_sol_prof_stable_zbot[lam], zarr_zbot_dict[lam] = sol_water[2], sol_water[6]


##plotting 
fig, ax = plt.subplots()
ax.set_ylabel('z[m]')
ax.set_xlabel('Eu')
ax.set_title('Eu Profile for Each Wavelength at zbot Where Ed = .001*E_d_0')
rgb_from_wavelengths = [wavelength_to_rgb.wavelength_to_rgb(k) for k in wavelengths]
colors = [(k[0]/255, k[1]/255, k[2]/255) for k in rgb_from_wavelengths] ##making a tuple from 0,1
for k,lam in enumerate(wavelengths) :
    ax.plot(Eu_sol_prof_stable_zbot[lam], zarr_zbot_dict[lam], label = wavelengths[k], color = colors[k])
ax.grid()
ax.legend(title='wavelengths[nm]')
    



























