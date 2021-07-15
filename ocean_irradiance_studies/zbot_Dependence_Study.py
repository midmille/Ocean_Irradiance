# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 05:52:04 2021

@author: Miles Miller
"""
"""
Zbot Dependence of Three Stream Irradiance Solution
----------------------------------------------------

--> Solves the three stream BVP for an array of points in ROMS output and for an
array of different zbots. 
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
## Importing other modules. 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import numpy as np




###########################zbot dependence sensitivity########################
#----------------------------------------------------------------------------#

##zbot study 
##Normalized version 
# wavelengths = [443]
# zbots =  np.arange(-15, -800,-5)
## a small call to see what the depths of the profiles from rep_points_index 
phy_profs = np.zeros((len(rep_points_index[:,0]), len(z_r[:,0,0])))
for k,(yi,xi) in enumerate(rep_points_index) : 
    phy_z_grid = z_r[:, int(yi), int(xi)]
    phy_profs[:,k] = phy_z_grid
diatom = roms_data.variables['diatom'][time_step_index,:,yi,xi] 
nanophyt = roms_data.variables['nanophytoplankton'][time_step_index,:,yi,xi]
phy_z_grid = z_r[:,yi,xi]
zbots_dict = {} ##fraction zbots of the stable zbot dict thus wavelength dep.
fract_zbot = np.arange(-.5,-5.25,-.1) ##to be multiplied by the zbot stable 
# fract_zbot = np.asarray([-.25, -.5, -1])
stable_zbot_dict = Make_Stable_zbot_Dict()
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
        solution = irr_phy_prof_shoot(lam, zbots_dict[lam][k], nanophyt, diatom, phy_z_grid, E_d_0, E_s_0, Nlayers=200,return_BVP=True)
        
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


def main():
    