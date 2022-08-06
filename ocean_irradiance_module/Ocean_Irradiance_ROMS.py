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
from multiprocessing import Pool
## can not import seapy at this point
# import os
## weird thing added because I am missing this value for some reason
# os.environ["PROJ_LIB"] = "C:/Users/Miles Miller/anaconda3/envs/pybasemap36/Library/share"
# import seapy 
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module import Ocean_Irradiance
from ocean_irradiance_module import Wavelength_To_RGB
import os
import sys
import pickle
import time
import matplotlib as mpl 
import matplotlib.pyplot as plt
if os.environ['CONDA_DEFAULT_ENV'] == 'ocean_irradiance': 
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature


def CI_alg(R_rs_b, R_rs_g, R_rs_r, lam_b, lam_g, lam_r): 
    """
    The CI algorithim from the NASA ocean color website.
   
    Parameters
    ----------
    R_rs_b : Float
        R_rs in the blue wavelength
    R_rs_g : Float 
        R_rs in the green wavelength
    R_rs_r : Float
        R_rs in the red wavelength
    lam_b : Float
        The wavelength in nm of blue. 
    lam_g : Float 
        The wavelength in nm of green. 
    lam_r : Float 
        The wavelength in nm of red.
    """

    return 

def OCx_alg(Rrs_443, Rrs_490, Rrs_510, Rrs_560, method='OC3M') :
    """
    Parameters
    ----------
    R_rs_b : Float
        R_rs in the blue wavelength
    R_rs_g : Float 
        R_rs in the green wavelength

    Returns
    -------
    chl_a : Float 
        Value of chl_a

    """
    OC3V = (0.2228, -2.4683, 1.5867, -0.4275, -0.7768) 
    OC4 = (0.3272, -2.9940, 2.7218, -1.2259, -0.5683)
    OC3M = (0.2424, -2.7423, 1.8017, 0.0015, -1.2280)
    OC4E = (0.3255, -2.7677, 2.4409, -1.1288, -0.4990)
    OC4O = (0.3325, -2.8278, 3.0939, -2.0017, -0.0257)
    

    ## the a values for OCx 0.2228	-2.4683	1.5867	-0.4275	-0.7768
    ## [These are the viirs values.]

    if method == 'OC3V': 
        a0, a1, a2, a3, a4 = OC3V
        R_rs_b = np.maximum(Rrs_443, Rrs_490) 
    elif method == 'OC4': 
        a0, a1, a2, a3, a4 = OC4
        R_rs_b = np.maximum(Rrs_443, Rrs_490) 
        R_rs_b = np.maximum(R_rs_b, Rrs_510) 
    elif method == 'OC3M': 
        a0, a1, a2, a3, a4 = OC3M
        R_rs_b = np.maximum(Rrs_443, Rrs_490) 
    elif method == 'OC4E': 
        a0, a1, a2, a3, a4 = OC4
        R_rs_b = np.maximum(Rrs_443, Rrs_490) 
        R_rs_b = np.maximum(R_rs_b, Rrs_510) 
    elif method == 'OC4O': 
        a0, a1, a2, a3, a4 = OC4O
        R_rs_b = np.maximum(Rrs_443, Rrs_490) 
        R_rs_b = np.maximum(R_rs_b, Rrs_510) 
    
    R_rs_g = Rrs_560
    
    
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
#    Q = 3
#    Q = 5
    ## This is Rrs sub surface. 
    ## See the following paper by Dutkiewicz et al. 2015.
    ## Modelling ocean -colour derived chlorophyll a 2018
    ## biogeosciences, 15, 613-630, 2018
    ## Rrs zero minus, subsurface
    R_rs_0m = R / Q
    
    ## HUGE AND IMPORTANT EDIT TO RRS
    ## Edit following Stephanie's paper for  
    R_rs = (0.52* R_rs_0m) / ( 1 - 1.7*R_rs_0m)
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


def Ocean_Irradiance_Field(PI, mask, ab_wat, ab_diat, ab_syn, chl_diatom, chl_nanophyt, 
                  z_r, coefficients, N=30, method='shootdown', zbot_arr = 'bathemetry'):
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
#    if method == 'shoot_down' or method == 'shoot_up' or method =='shoot_fp' or method == 'scipy':
#        irr_arr = np.zeros((N,nyi,nxi,4)) * np.nan 
#    if method == 'Dut':
#        irr_arr = np.zeros((N,nyi,nxi,4)) * np.nan 

    irr_arr = np.zeros((N,nyi,nxi,4)) * np.nan 

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
                if zbot_arr == 'bathemetry':
                    hbot = z_r0[0]
                else: 
                    hbot = zbot_arr[j,i]

                assert len(chl_diatom) == len(z_r0)
                
                phy_profs = np.zeros((len(z_r0),2))
                phy_profs[:,0] = chl_diatom[:,j,i]
                phy_profs[:,1] = chl_nanophyt[:,j,i]
                #print(phy_profs)
                
                a = np.array([ab_diat[0], ab_syn[0]])
                b = np.array([ab_diat[1], ab_syn[1]])
                
                phy = Ocean_Irradiance.Phy(z_r0, phy_profs, [esd('Diat'), esd('Syn')], a, b)    

                ocean_irr_sol = Ocean_Irradiance.ocean_irradiance(PI, 
                                                                  hbot, 
                                                                  ab_wat, 
                                                                  method = method, 
                                                                  phy=phy, 
                                                                  N=N)

                    
                ## Ed
                irr_arr[:,j,i,0] = ocean_irr_sol[0]
                ## Es
                irr_arr[:,j,i,1] = ocean_irr_sol[1]
                ## Eu
                irr_arr[:,j,i,2] = ocean_irr_sol[2]
                ## z_irr
                irr_arr[:,j,i,3] = ocean_irr_sol[3]

                
    return irr_arr


def Irradiance_Run(R_nc, PI, nstp, N, savefile, method, wavelengths):
    
    
    ## Python calculated Eu at surface 
    ##---------------------------------
    ## Checking if save file exists
    ## if it doesn't exist then redo calculation
    if os.path.exists(savefile) == False:
        ## User input required so not to overwrite or redo unnecessarily... 
        #y_or_n = input('File does not exist, continue with calculation? [y/n] ')
        #if y_or_n == 'n': 
        #    sys.exit('Stopping...')
        #elif y_or_n == 'y':
        #    print('ok, statrting calculations...')
            
        #mask = np.ones((R_nc.nyi, R_nc.nxi))
        irr_field_py = {}
        for lam in wavelengths:
            print('Current Wavelength:', lam)
            irr_field_py[lam] = Ocean_Irradiance_Field(
                                              PI,
                                              R_nc.maskr, 
                                              abscat(lam, 'water'), 
                                              abscat(lam, 'Diat'), 
                                              abscat(lam, 'Syn'), 
                                              R_nc.chl_diatom[nstp,:,:,:], 
                                              R_nc.chl_nanophyt[nstp,:,:,:], 
                                              R_nc.z_r[nstp,:,:,:], 
                                              PI.coefficients,
                                              N= N, 
                                              method = method)
        
        pickle.dump(irr_field_py, open(savefile, "wb"))
        print('Python calculation complete and saved')
        
    ## if the save file does exist then just load it. gity 
    elif os.path.exists(savefile) == True:
        print(f'Irradiance save file exists! Loading python calculated irradiance field from file "{savefile}"...')
        irr_field_py = pickle.load(open(savefile,'rb'))
        print('Yay, file loaded :)')

    return irr_field_py


def run_irradiance_one(args): 
    lam, PI, R_nc, N, method, nstp = args
    irr_field_py = Ocean_Irradiance_Field(
                                            PI,
                                            R_nc.maskr, 
                                            abscat(lam, 'water'), 
                                            abscat(lam, 'Diat'), 
                                            abscat(lam, 'Syn'), 
                                            R_nc.chl_diatom[nstp,:,:,:], 
                                            R_nc.chl_nanophyt[nstp,:,:,:], 
                                            R_nc.z_r[nstp,:,:,:], 
                                            PI.coefficients,
                                            N= N, 
                                            method = method)
        
        
    return lam, irr_field_py 
 

def Irradiance_Run_Parallel(R_nc, PI, nstp, N, savefile, method, wavelengths):

   
    ## [The number of processors is equal to the number of wavelengths.]
    num_procs = len(wavelengths)
    
    ## Python calculated Eu at surface 
    ##---------------------------------
    ## Checking if save file exists
    ## if it doesn't exist then redo calculation
    if os.path.exists(savefile) == False:
        irr_field_py = {}

        ## [Run in parallel.]
        args=[]
        for lam in wavelengths: 
            args.append([lam, PI, R_nc, N, method, nstp])

        process_pool = Pool(num_procs)
        result = process_pool.map(run_irradiance_one, args)
        for res in result: 
            irr_field_py[res[0]] = res[1]
        process_pool.close()

       
        pickle.dump(irr_field_py, open(savefile, "wb"))
        print('Python calculation complete and saved')
        
    ## if the save file does exist then just load it. gity 
    elif os.path.exists(savefile) == True:
        print("The file already exists! To load file please do not run with parallel flag.")

    return 



def Plot_Irradiance_Field(R_nc, irr_field, nstp, Ed0, Es0):
    """
    This function is to plot the irradiancee field and other values from an irradiance
    run over a ROMS out put. The things it plots are as follows:
    1) Log of Nanophytoplankton Concentration.
    2) Log of Diatom Concentrations. 
    3) The Rrs values from irradiance model for three given wavelengths.
    4) The OCx Chl-a calculated from the Rrs.
    """

    ## Getting lat, lon, nano, and diat from R_nc.
    ## The negative one is to get the surface value.
    chla_diat = R_nc.chl_diatom[nstp,-1,:,:]
    chla_nano = R_nc.chl_nanophyt[nstp,-1,:,:]
    lat = R_nc.lat_rho
    lon = R_nc.lon_rho
    
    ## The figure projection of lat/lon.
    proj = ccrs.PlateCarree()
    ## makes the color bar on plots a bit smaller if necessary.
    cbar_shrink = 1 
    
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2,8)
    ax00 = fig.add_subplot(gs[0, 0:2], projection=proj)
    ax01 = fig.add_subplot(gs[0, 2:4], projection=proj)
    ax02 = fig.add_subplot(gs[0, 4:6], projection=proj)
    ax10 = fig.add_subplot(gs[1, 0:2], projection=proj)
    ax11 = fig.add_subplot(gs[1, 2:4], projection=proj)
    ax12 = fig.add_subplot(gs[1, 4:6], projection=proj)
    ## This ax is for the rainbow legend plot.
    ax_tall = fig.add_subplot(gs[0:2, 7])

    ## The Rrs plots list. 
    Rrs_axs = [ax10, ax11, ax12]

    for k,lam in enumerate(wavelengths):
         
        ## Getting Eu at the surface
        Eu_surf = irr_field[lam][-1, :,:, 2]
        ## Calculating Rrs
        Rrs = R_RS(Ed0, Es0, Eu_surf) 
        if lam == 443:
            Rrs_b = Rrs
        if lam == 551:
            Rrs_g = Rrs
        
        ## Plotting Rrs map
        ax = Rrs_axs[k]
        ## Plotting the land and other features to make pretty.
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
        #ax.gridlines()
        im = ax.pcolormesh(lon, lat, Rrs, cmap='nipy_spectral', 
                           transform=ccrs.PlateCarree())#, vmax=vmax, vmin=vmin )
        fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = r'Rrs [$\mathrm{sr}^{-1}$]')
        ax.set_title(r'Rrs $\lambda=$ ' + f'{lam}')  
        ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
        ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))
        
        ### Including the absorbtion/scat vals
        #ab_wat = abscat(lam,'water')
        #ab_diat = abscat(lam, 'Diat')
        #ab_syn = abscat(lam, 'Syn')
        ### Starting points for text.
        #x = lon[int(lon.shape[1] * (2/3))]
        #y = lat[int(lat.shape[0] * (2/3))]
        ### Text sep is 1/10 of the total lat
        #text_sep = lat[int(lat.shape[0] * (1/10))]
        #font_size = 10
        #ax.text(x,y,r'$a_{wat} = $ %.4f'%ab_wat[0], fontsize = font_size)
        #y = y - text_sep
        #ax.text(x,y, r'$b_{wat} = $ %.4f'%ab_wat[1], fontsize = font_size)
        #y = y - text_sep
        #ax.text(x,y, r'$a_{diat} = $ %.4f'%ab_diat[0], fontsize = font_size)
        #y = y - text_sep
        #ax.text(x,y, r'$b_{diat} = $ %.4f'%ab_diat[1], fontsize = font_size)
        #y = y - text_sep
        #ax.text(x,y, r'$a_{syn} = $ %.4f'%ab_syn[0], fontsize = font_size)
        #y = y - text_sep
        #ax.text(x,y, r'$b_{syn} = $ %.4f'%ab_syn[1], fontsize = font_size)
    
    ## The OCx calculation
    chla_ocx = OCx_alg(Rrs_b, Rrs_g)
    
    ## list of different chla plots.
    chla_axs = [ax00, ax01, ax02]
    chla_list = [chla_ocx, chla_diat, chla_nano] 
    chla_labels = ['OCx Chl-a from Irradiance Rrs', 'NEMURO Diatoms', 'NEMURO Nanophytoplankton']
    for k, ax in enumerate(chla_axs):
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
        #ax.gridlines()
        im = ax.pcolormesh(lon, lat, chla_list[k], cmap='nipy_spectral', 
                           transform=ccrs.PlateCarree())#, vmax=vmax, vmin=vmin )
        fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = r'Chl-a [mg Chl-a $\mathrm{m}^{-3}$')
        ax.set_title(chla_labels[k])  
        ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
        ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))

    ## Plotting the rainbow with wavelength indicators.
    ## number of wavelength points.
    ## The R stands for rainbow to differntiate between lams used for irradiance. 
    R_Nlam = 750 - 380
    R_lams = np.linspace(381, 750, R_Nlam)
    R_lams_rgb = np.zeros((R_Nlam, 3))
    for k, R_lam in enumerate(R_lams):
        R_lams_rgb[k] = Wavelength_To_RGB.wavelength_to_rgb(R_lam)
    len_x = 100 
    ## The rgb matrix
    R_rgb = np.zeros((R_Nlam, len_x, 3)) 
    for k in range(len_x):
        R_rgb[:,k,:] = R_lams_rgb
    ## The actual plotting of the wavelength rainbow.
    ax_tall.imshow(R_rgb, origin='lower', aspect=('auto'))
    ## The indexes of the irr wavelengths.
    indexes = []
    for lam in wavelengths:
        index = lam - 380
        ax_tall.plot([0, len_x-1], [index, index], 'k') 
        indexes.append(index)
    ax_tall.set_yticks(indexes)
    ax_tall.set_yticklabels(wavelengths)
    ax_tall.set_xticks([])

    ##plotting some date text. The date from ocean_time. 
    #plt.figtext(.02,.95,'ROMS Date and Time: \n {}'.format(roms_date[0]))
    ##plotting some units text for absorbtion and scattering
    #plt.figtext(.02,.04, 'a_wat [m^-1] \nb_wat [m^-1] \na_phy [m^2 / mg Chl] \nb_phy [m^2 / mg Chl]')

    fig.show()
    ###The rainbow with markers of wavelengths used 
    #plt.subplot(3,10,(10,30)) ##making it a small top thing at the top of the array
    #NP = 750 - 380
    #P = np.linspace(381,750,NP)
    #P_rgb = np.zeros((NP,3))
    #for k in range(NP): 
        #P_rgb[k] = wavelength_to_rgb.wavelength_to_rgb(P[k])
   # 
    #len_x = 100
    #len_y = len(P)
    #T = np.zeros((len_y,len_x,3)) ##slightly confusing but matrix is oriented this way 
    #for i in range(len_x):
        #T[:,i,:] = P_rgb / 255
   # 
    #plt.imshow(T,origin = 'lower', aspect=('auto'))
    ## plt.axis('off')
    #indexes = []
    #for lam in wavelengths: 
        #index = lam - 380
        #plt.plot([0,len_x-1], [index,index], 'k')
        ## plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
        #indexes.append(index)
    #plt.yticks(indexes, wavelengths)
    #plt.xticks([])
    
    ##Plotting the composite rgb_grid image 
    #plt.subplot(3,5,2)
    #plt.imshow(rgb_grid, origin=('lower'), aspect=('auto'))
    #plt.title('True Color \n Wavelength Composite')
    #plt.axis('off')
    # plt.colorbar()
    
    #plt.pcolormesh(diatom, norm = mpl.colors.LogNorm())
    #plt.title('Log Plot Diatom Surface Concentration \n [{}]'.format(diat_units))
    #plt.colorbar()
    #plt.axis('off')
   # 
    #plt.subplot(3,5,4)
    #plt.pcolormesh(nanophyt, norm = mpl.colors.LogNorm())
    #plt.title('Log Plot Nanophytoplankton(Syn) \n Surface Concentration \n [{}]'.format(nano_units))
    #plt.colorbar()
    #plt.axis('off')
    
    
    return


    return
    
if __name__ == '__main__':
    
    import argparse
    ## User Modules
    from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    parser = argparse.ArgumentParser(description='Ocean Irradiance ROMS Wrapper')
    parser.add_argument('romsnc_file', help = "Complete Path to ROMS nc file" )
    parser.add_argument('savefile', help = "The name of the pickle file in which the result will be stored" )
    parser.add_argument('method', help = "methods are: shootdown, shootup, scipy, dutkiewicz" )
    # parser.add_argument('dest_file', help='Path to Destination Directory. Saved as pickle')
    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
    parser.add_argument('--parallel', action='store_true', help="Runs the wavelengths on seperate processors")
    args = parser.parse_args()
    

    PI = Param_Init()

    PI.pt1_perc_zbot = True
    PI.pt1_perc_phy = True

    PI.grid = 'log'

    R_nc = ROMS_netcdf(args.romsnc_file)

    ## updating default to wavelengths necessary for the OCx_algorithim
    wavelengths = PI.wavelengths

    ## settting time step index
    time_step_index = 1
    
    N = 30

    if args.parallel: 
        t1 = time.perf_counter()
        Irradiance_Run_Parallel(R_nc, PI, time_step_index, N, args.savefile, args.method, wavelengths)
        t2 = time.perf_counter()
        print('Compute time: ', t2-t1)
    else: 
        t1 = time.perf_counter()
        irr_field = Irradiance_Run(R_nc, PI, time_step_index, N, args.savefile, args.method, wavelengths)
        if args.plot:
            Plot_Irradiance_Field(R_nc, irr_field, time_step_index, PI.Ed0, PI.Es0)
        t2 = time.perf_counter()
        print('Compute time: ', t2-t1)
        

   
