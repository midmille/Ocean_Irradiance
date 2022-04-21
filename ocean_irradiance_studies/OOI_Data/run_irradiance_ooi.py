"""
Created: March 9, 2022 7:48am
Author: Miles D. Miller, University of California Santa Cruz

This file is for the purpose of running the irradaince light model using 
OOI data as well as comparing the model output to in situ data.

To run: 
    --> Be sure to be in the ooi environment. 
        -https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python
"""

## [Import Statements]
## [User Modules]
import ooi_data_functions as ODF 
from ocean_irradiance_module import Ocean_Irradiance as OI
import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module import Wavelength_To_RGB

## [External Modules]
import numpy as np
import pickle
from ooi_data_explorations.data_request import data_request
import matplotlib.pyplot as plt
import scipy


def Run_Irradiance(N, wavelengths, phy_type, depth_profs, dt_profs, chla_profs, opt_bs_profs, optaa_dat):
    """
    This function runs the irradiance model over all the different chla profiles.
    
    It runs the irradiance model twice, the first time using the coefficients provided by Dutkiewicz using OOI
    chla profiles, and the second time using OOI absorption and back scatter profiles.

    Parameters
    ----------

    Returns 
    -------

    """
    ## [The irradiance parameters object.]
    PI = Param_Init()

    ## [The number of profiles.]
    #N_profs = len(depth_profs)
    N_profs = 1

    ## [Irradiance field dictionary.]
    irr_field = {}
    irr_field_ab = {}

    ## [Loop over the wavelengths.]
    for lam in wavelengths:
        ## [The irradiance array.]
        irr_arr = np.zeros((N, N_profs, 4))
        irr_arr_ab = np.zeros((N, N_profs, 4))

        ## [Loop over the OOI profiles.]
        for k in range(N_profs):
            ## [Make the depth negative.]
            depth = -depth_profs[k]
            ## [Make the depth profile start at zero.]
            ## [Create the phytoplankton object.]
            phy = OI.Phy(depth, chla_profs[k], esd(phy_type), 
                         abscat(lam, phy_type, C2chla='default')[0], 
                         abscat(lam, phy_type, C2chla='default')[1])
            ocean_irr_sol = OIS.ocean_irradiance_shubha(depth[0], 
                                                        PI.Ed0+PI.Es0, 
                                                        abscat(lam, 'water'), 
                                                        PI.coefficients, 
                                                        phy=phy, 
                                                        CDOM=None, 
                                                        N=N, 
                                                        pt1_perc_zbot=False, 
                                                        pt1_perc_phy=False)
            ## [Storing output into array.]
            irr_arr[:,k,0] = ocean_irr_sol[0]
            irr_arr[:,k,2] = ocean_irr_sol[1]
            irr_arr[:,k,3] = ocean_irr_sol[2] 
            a_irr = ocean_irr_sol[3]
            b_b_irr = ocean_irr_sol[4]

            ## [Now using the OOI absorption and scattering.]
            ## [The optaa data.]
            optaa_depth_dat = optaa_dat.variables['depth'].data
            optaa_dt_dat = optaa_dat.variables['time'].data
            optaa_abs_dat = optaa_dat.variables['optical_absorption'].data
            optaa_wavelength_dat = optaa_dat.variables['wavelength_a'].data

            dt_lbnd = dt_profs[k].data[0] 
            dt_ubnd = dt_profs[k].data[-1]
            prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, optaa_dt_dat)

            optaa_dt_prof = optaa_dt_dat[prof_mask]
            optaa_depth_prof = optaa_depth_dat[prof_mask]
            optaa_abs_prof = optaa_abs_dat[prof_mask]
            optaa_wavelengths = optaa_wavelength_dat[prof_mask]

            optaa_depth = -np.squeeze(optaa_depth_prof)
            ## [Must get the wavelength index.]
            lam_i = ODF.Get_Wavelength_Index(optaa_wavelengths[0,:], lam)

            ## [The squeese makes the array 1-D.]
            a = np.squeeze(optaa_abs_prof[:, lam_i])
            b = opt_bs_profs[k]

            ## [Now running the irr model using abs and scat from ooi.]
            ocean_irr_sol_ab = OIS.ocean_irradiance_shubha_ab(depth[0], 
                                                              PI.Ed0+PI.Es0, 
                                                              PI.coefficients, 
                                                              optaa_depth.data,  ## [The z_a is optta depth.]
                                                              a.data, 
                                                              depth.data, ## [The z_b is the flort depth.]
                                                              b.data,
                                                              CDOM=None, 
                                                              N=N, 
                                                              pt1_perc_zbot=False)


            ## [Storing output into array.]
            irr_arr_ab[:,k,0] = ocean_irr_sol_ab[0]
            irr_arr_ab[:,k,2] = ocean_irr_sol_ab[1]
            irr_arr_ab[:,k,3] = ocean_irr_sol_ab[2] 
        
        ## [Store array to dictionary.]
        irr_field[lam] = irr_arr
        irr_field_ab[lam] = irr_arr_ab

    return irr_field, irr_field_ab

def OOI_Abs_Scat(prof_index, lam, depth_profs, dt_profs, opt_bs_profs, optaa_dat):
    """
    This function gets the OOI absorption and scattering for a given profile index and data sets.
    """

    k = prof_index

    ## [The optaa data.]
    optaa_depth_dat, optaa_dt_dat, optaa_abs_dat, optaa_wavelength_dat = ODF.Get_Optaa_Dat(optaa_dat)

    dt_lbnd = dt_profs[k].data[0] 
    dt_ubnd = dt_profs[k].data[-1]
    prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, optaa_dt_dat)

    optaa_dt_prof = optaa_dt_dat.data[prof_mask]
    optaa_depth_prof = optaa_depth_dat.data[prof_mask]
    optaa_abs_prof = optaa_abs_dat.data[prof_mask]
    optaa_wavelengths = optaa_wavelength_dat.data[prof_mask]

    optaa_depth = np.squeeze(optaa_depth_prof)
    ## [Must get the wavelength index.]
    lam_i = ODF.Get_Wavelength_Index(optaa_wavelengths[0,:], lam)
    lam_ooi = optaa_wavelengths[0,lam_i]

    ##[Labeling the ooi z-grids and a, b_b.]
    z_a = optaa_depth
    ## [The squeese makes the array 1-D.]
    a = np.squeeze(optaa_abs_prof[:, lam_i])
    z_b_b = depth_profs[k].data
    b_b = opt_bs_profs[k].data

    return z_a, a, z_b_b, b_b


def Irr_OOI_Abs_Scat(PI, N, lam, prof_idex, phy_type, depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs, optaa_dat): 
    """
    This function compares the absorption and backscatter as calculated by the irradiance model with Dutkiewicz
    coefficients to the values given by OO.
    """

    
    phy = OI.Phy(depth, chla_profs[k], esd(phy_type), 
                 abscat(lam, phy_type, C2chla='default')[0], 
                 abscat(lam, phy_type, C2chla='default')[1])

    CDOM2C = 0
    CDOM_dens = OI.CDOM_dens(depth, cdom_profs[k], CDOM2C, lam)
            
    z_irr, a_irr, b_irr, b_b_irr  = OIS.ocean_irradiance_two_stream_ab(depth[0], 
                                                                   abscat(lam, 'water'), 
                                                                   N,
                                                                   phy=phy, 
                                                                   CDOM_dens=CDOM_dens, 
                                                                   pt1_perc_zbot=False, 
                                                                   pt1_perc_phy=False)

    ## [Get the ooi absorption and scattering, along with their corresponding grids.]
    z_a_ooi, a_ooi, z_b_b_ooi, b_b_ooi = OOI_Abs_Scat(prof_index, lam, depth_profs, dt_profs, opt_bs_profs, optaa_dat) 


    return  z_irr, a_irr, b_b_irr, z_a_ooi, a_ooi, z_b_b_ooi, b_b_ooi, lam_ooi


def Abs_Est_Species(PI, N, lam, prof_index, phy_species, depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs, optaa_dat): 
    """
    This function is for the estimation of the ratio of different species that could compose the
    resulting observed absorption profile from OOI. The total observed phytoplankton concentration will be
    taken from the FLORT instrument and the observed absorption will be taken from the OPTAA instrument. 

    This estimation will be performed as a least squares problem where the different points in the z-axis will
    be considered the data points sucht that the system is overdetermined. The downside to this is that is assumes the ratios
    of different phytoplankton will be constant within a single profile. Perhaps down the line the water column can be split up into different
    sections and we can solve for the phytoplankton in each different overdetermined section. 

    The matrix problem is as follows: 

    Ax = y, where 

    ## [The effective absorption coefficient for each species.]
    A = [[ a_phy1, a_phy2, a_phy3, a_phy4, ....], 
         [ a_phy1, a_phy2, a_phy3, a_phy4, ....],
         [    :  ,   :   ,  :    ,   :   , ....], 
         [    :  ,   :   ,  :    ,   :   , ....]]

    ## [The ratios of the different phytoplankton species, this is what we are solving for.]
    x = [ r_phy1, r_phy2, r_phy3, r_phy4, ....]^T 

    ## [The effective absorption given by observed OOI absorption divided by the 
        observed OOI chla and then subtracted by the dutkiewicz a_wat.]
    y = [a_ooi(z1), a_ooi(z2), a_ooi(z3), a_ooi(z4), a_ooi(z4) ... ]^T

    Note: The OOI absorption comes from the OPTAA instrument while the chla comes from the
    FLORT instrument. Thus it is necessary to interpolate them to the same grid. 

    Parameters
    ----------

    Returns 
    -------

    """

    ## [The number of species.]
    N_phy = len(phy_species)

    ## [Get the ooi absorption and scattering, along with their corresponding grids.]
    z_a_ooi, a_ooi, z_bb_ooi, bb_ooi = OOI_Abs_Scat(prof_index, lam, depth_profs, dt_profs, opt_bs_profs, optaa_dat) 
    
    ## [Get the arrays for the given prof_index.]
    ## [This is the chla grid and the grid to be used for everything.]
    z_chla = depth_profs[prof_index].data
    chla = chla_profs[prof_index].data
    cdom = cdom_profs[prof_index].data

    ## [The chla and b_b ooi should be on same grid.]
    ## [Interpolate the ooi chla grid, and ooi b_b to the ooi absorption grid.]
    z = z_a_ooi
    chla = np.interp(z, z_chla, chla)
    bb_ooi = np.interp(z, z_bb_ooi, bb_ooi)

    ## [The number of vertical levels in the ooi grid.]
    N = len(z)

    ## [Construct the empty matrix system.]
    A = np.zeros((N,N_phy))
    b = np.zeros(N)

    ## [Loop over the phytoplankton species.]
    for k in range(N_phy): 
        ## [Each column of a is the single effective abs/scat coefficient of a species, 
        ##  each column corresponds to a different species.]
        A[:,k] = (abscat(lam, phy_species[k], C2chla='default')[1] * chla)

    ## [The rhs b, is simply the a_ooi divided by chla minus the a_wat from Dut.]
    b = bb_ooi - abscat(lam, 'water', C2chla='default')[1]
        
    ## [Solving the least square system for the phy ratios.]
    #x = np.linalg.lstsq(A, b)
    x = scipy.optimize.nnls(A,b)

    #phy = OI.Phy(depth, chla_profs[k], esd(phy_type), 
    #             abscat(lam, phy_type, C2chla='default')[0], 
    #             abscat(lam, phy_type, C2chla='default')[1])

    #CDOM2C = 0
    #CDOM_dens = OI.CDOM_dens(depth, cdom_profs[k], CDOM2C, lam)

    ## [Plot the resulting absorption from the calculated ratios of phytoplankton.]
    fig, ax = plt.subplots()

    ## [Plot the 
    a_fit = A[0,:]@x[0]
    for k in range(N_phy): 
        a_f = a_fit[k] * chla
    ax.plot(


    return A, b, x 


def Plot_Irr_OOI_Abs_Scat(PI, prof_index, wavelengths, N, phy_types, depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs, optaa_dat): 
    """

    """

    ## [Get the colors such that they match the wavelengths.] 
    colors = [Wavelength_To_RGB.wavelength_to_rgb(wavelength) for wavelength in wavelengths]
    
    N_lam = len(wavelengths)
    figa, axas = plt.subplots(ncols=N_lam, nrows=1)
    figb, axbs = plt.subplots(ncols=N_lam, nrows=1)
    for k, lam in enumerate(wavelengths):
        ## [The different plots correspond to a given wavlength]
        axa = axas[k]
        axb = axbs[k]
        for i, phy_type in enumerate(phy_types):
            ## [Getting the absorption and scattering values.]
            z_irr, a_irr, b_b_irr, z_a_ooi, a_ooi, z_b_b_ooi, b_b_ooi, lam_ooi = Irr_OOI_Abs_Scat(PI,
                                                                                                 N,
                                                                                                 lam,
                                                                                                 prof_index,
                                                                                                 phy_type,
                                                                                                 depth_profs, 
                                                                                                 dt_profs,
                                                                                                 chla_profs, 
                                                                                                 cdom_profs,
                                                                                                 opt_bs_profs, 
                                                                                                 optaa_dat)
            
            ## [Plotting Absorption comparison.]
            ## [Only plot the ooi abs, scat for i==0, since they are indp. of wavelength.]
            if i == 0: 
                axa.plot(a_ooi, z_a_ooi, ':', label=f'Abs OOI OPTAA {lam_ooi}')
            axa.plot(a_irr, z_irr, label=f'Abs Irr {phy_type}')
            axa.set_title(f'Absorption {lam}')
            axa.set_ylabel('Z[m]')
            axa.set_xlabel('Abs [m^-1]')
            axa.legend()
            axa.grid()
    
            ## [Plotting the scattering comp.]
            ## [Only plot the ooi abs, scat for i==0, since they are indp. of wavelength.]
            if i == 0: 
                axb.plot(b_b_ooi, z_b_b_ooi, ':', label=f'Scat OOI FLORT {lam_ooi}')
            axb.plot(b_b_irr, z_irr, label=f'Scat Irr {phy_type}')
            axb.set_title(f'Back Scattering {lam}')
            axb.set_ylabel('Z[m]')
            axb.set_xlabel('Scat [m^-1]')
            axb.legend()
            axb.grid()
 
    
    figa.show()
    figb.show()

    return 


def Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_wavelengths, spkir_wavelength_index,irr_field, irr_field_ab, spkir_dat, site, assembly, method): 

    """
    This function plots downwelling irradiance stream solved for using the chla profiles
    from ooi. Then it plots the spkir data that corresponds to the solved irradiance profiles.

    Parameters
    ----------
    prof_index: Integer
        The index of the profile to be plotted. 
    wavelengths: List, Integer
        The list of wavelengths, the wavelengths of the ooi SPKIR might differ from the
        wavelengths used for the coefficients in the irradiance model. Might has to pass
        the spkir wavelengths eventually. 
    irr_field: Dictionary, 4-D Array
        The irradiance field dictionary with wavelengths as keys and the 4-D irradiance array
        as the values.
    spkir_depth_profs: List, 1-D Array
        The list of spkir depth profiles. The depth profiles are positive with the bottom most
        depth at index 0 and the surface at index -1. 
    spkir_dt_profs: List, 1-D Array
        The list of the times associated with each depth and spkir data point. 
    spkir_profs: List, 2-D Array
        The list of spkir profiles. The first index of each lsit value is the profile and the
        second index corresponds to the wavelenth.
    """

    ## [The number of wavelengths.]    
    N_lam = len(wavelengths)

    ## [spkir depth stuff.]
    spkir_depth_dat = spkir_dat.variables['depth'].data
    spkir_dt_dat = spkir_dat.variables['time'].data
    spkir_dat = spkir_dat.variables['spkir_abj_cspp_downwelling_vector'].data

    dt_lbnd = dt_profs[prof_index].data[0] 
    dt_ubnd = dt_profs[prof_index].data[-1]
    prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, spkir_dt_dat)

    spkir_dt_prof = spkir_dt_dat[prof_mask]
    spkir_depth_prof = spkir_depth_dat[prof_mask]
    spkir = spkir_dat[prof_mask]

    depth = -spkir_depth_prof 

    fig, ax = plt.subplots()
    ## [Get the colors such that they match the wavelengths.] 
    colors = [Wavelength_To_RGB.wavelength_to_rgb(wavelength) for wavelength in wavelengths]
    
    ## [Loop over the irradiance wavelengths.]
    for k, lam in enumerate(wavelengths): 
        irr_arr = irr_field[lam]
        irr_arr_ab = irr_field_ab[lam]
        ## [ Plotting the downward direct profile. 0 is downward irradiance, 3 is depth.]
        ## [Also must multiply by surface value of spkir prof, since irr_arr is normalized.]
        ax.plot(spkir[:,spkir_wavelength_index[k]].data[-1] * irr_arr[:, prof_index, 0], irr_arr[:, prof_index, 3], ':', label=f'Irr {lam}', color=colors[k])
        ax.plot(spkir[:,spkir_wavelength_index[k]].data[-1] * irr_arr_ab[:, prof_index, 0], irr_arr_ab[:, prof_index, 3], '-', label=f'Irr ooi ab {lam}', color=colors[k])

    for k,i in enumerate(spkir_wavelength_index):
        ## [Plotting the spkir profile.]
        ## [Make the 1m avg grid.]
        depth_avg = np.arange(depth[0], depth[-1], 1)
        spkir_avg = ODF.Grid_Average_Profile(depth, spkir[:,i], depth_avg)
        #ax.plot(spkir[:, i], depth, '--', label=f'OOI SPKIR {lam}', color=colors[k])
        ax.plot(spkir_avg, depth_avg, '--', label=f'OOI SPKIR {lam}', color=colors[k])

    ## [Labels.]
    #ax.set_ylabel(f"Z [{depth_dat.attrs['units']}]")
    ax.set_ylabel(f"Z [m]")
    #ax.set_xlabel(f"Downwelling Spectral Irradiance {spkir_dat.attrs['units']}")
    ax.set_xlabel(f"Downwelling Spectral Irradiance")
    #ax.set_title(f"OOI SPKIR Profile and Irradiance Model \n OOI SPKIR Profile Date: {date_time[0]} to {date_time[-1]}")
    ax.set_title(f"OOI SPKIR Profile and Irradiance Model")
    ## [Putting some identifying text on the figure.]
    ## [10% up the vertical location]
    txt_y = ax.get_ylim()[1] + 0.5 * ax.get_ylim()[0] 
    ## [10% of the horizontal location.]
    txt_x = ax.get_xlim()[0] + 0.2 * ax.get_xlim()[1]
    ## [The change in txt location in vertical.]
    txt_dz = 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ## [Adding the txt.]
    ax.text(txt_x, txt_y, f'SITE: {site}')   
    ax.text(txt_x, txt_y+txt_dz, f'ASSEMBLY: {assembly}')   
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: SPKIR')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.legend(title='Wavelengths [nm]')
    ax.grid()

    fig.show()

    return 

def Plot_OOI_Abs_Wavelength_Time(optaa_dat, depth_ref, start, stop, site_name, assembly, method):
    """
    This plot creates a plot similiar to the one Chris Wingard sent in his email on 4/9/2022. 
    """

    ## [The optaa data.]
    depth_profs, dt_profs, abs_profs, wavelength_profs = Get_Optaa_Profiles(optaa_dat, abs_cor=True)

    ## [The number of profiles.]
    N = len(depth_profs)
    ## [The number of wavelengths.]
    optaa_wavelengths = np.squeeze(wavelength_profs[0][0,:])
    N_lam = len(optaa_wavelengths)

    ## [The absorption in time.]
    abs_t = np.zeros((N, N_lam))
    ## [Loop over the profiles.]
    for k in range(N): 
        depth = depth_profs[k].data
        ## [The depth index.]
        d_i = np.argmin(abs(depth-depth_ref)) 
        ## [The absorption in time.]
        abs_t[k,:] = abs_profs[k][d_i, :] 

    ## [Plotting.]
    fig, ax = plt.subplots()
    
    ## [The actual plotting]
    ax.plot(optaa_wavelengths, np.transpose(abs_t), 'b', linewidth=.7)  

    ## [Labels.]
    ax.set_xlabel(f"Wavelength [{optaa_dat.variables['wavelength_a'].attrs['units']}]")
    ax.set_ylabel(f"Absorption [{optaa_dat.variables['optical_absorption'].attrs['units']}]")
    ax.set_title(f"OOI OPTAA optical_absorption Multiple Profiles Single Depth ({depth_ref} [m]) \n START: {start} to STOP: {stop}")
    ## [Putting some identifying text on the figure.]
    ## [10% up the vertical location]
    txt_y = ax.get_ylim()[0] + 0.5 * (ax.get_ylim()[1]  - ax.get_ylim()[0])
    ## [10% of the horizontal location.]
    txt_x = ax.get_xlim()[0] + 0.6 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    ## [The change in txt location in vertical.]
    txt_dz = 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ## [Adding the txt.]
    ax.text(txt_x, txt_y, f'SITE: {site_name}')   
    ax.text(txt_x, txt_y+txt_dz, f'ASSEMBLY: {assembly}')   
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: OPTAA')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.grid()
    
    fig.show()

    return
    

def Plot_OOI_Abs_Wavelength_Prof(optaa_dat, prof_index, start, stop, site_name, assembly, method):
    """
    This plot creates a plot similiar to the one Chris Wingard sent in his email on 4/9/2022. 
    The main difference between this function and Plot_OOI_Abs_Wavelength_Time() is that this function
    plots the absorption against wavelength for all depths at single wavelengths. 
    """

    ## [The optaa data.]
    depth_profs, dt_profs, abs_profs, wavelength_profs = Get_Optaa_Profiles(optaa_dat, abs_cor=True)

    ## [The number of profiles.]
    N = len(depth_profs[prof_index])
    ## [The number of wavelengths.]
    optaa_wavelengths = np.squeeze(wavelength_profs[0][0,:])
    N_lam = len(optaa_wavelengths)

    ## [The absorption in profile.]
    abs_p = np.zeros((N, N_lam))

    ## [Loop over the profiles.]
    for k in range(N): 
        ## [The absorption in the profile]
        abs_p[k,:] = abs_profs[prof_index][k, :] 

    print( 'Number of lines:', N)

    ## [Plotting.]
    fig, ax = plt.subplots()

    ## [The plotting.]
    ax.plot(optaa_wavelengths, np.transpose(abs_p), 'b', linewidth=.7)  

    ## [Labels.]
    ax.set_xlabel(f"Wavelength [{optaa_dat.variables['wavelength_a'].attrs['units']}]")
    ax.set_ylabel(f"Absorption [{optaa_dat.variables['optical_absorption'].attrs['units']}]")
    #ax.set_title(f"OOI SPKIR Profile and Irradiance Model \n OOI SPKIR Profile Date: {date_time[0]} to {date_time[-1]}")
    dt_sec = dt_profs[prof_index].data.astype('datetime64[s]')
    ax.set_title(f"OOI OPTAA optical_absorption Single Profile All Depths \n START: {dt_sec[0]} to STOP: {dt_sec[-1]}")
    ## [Putting some identifying text on the figure.]
    ## [10% up the vertical location]
    txt_y = ax.get_ylim()[0] + 0.7 * (ax.get_ylim()[1]  - ax.get_ylim()[0])
    ## [10% of the horizontal location.]
    txt_x = ax.get_xlim()[0] + 0.6 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    ## [The change in txt location in vertical.]
    txt_dz = 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ## [Adding the txt.]
    ax.text(txt_x, txt_y, f'SITE: {site_name}')   
    ax.text(txt_x, txt_y+txt_dz, f'ASSEMBLY: {assembly}')   
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: OPTAA')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.grid()
    
    fig.show()

    return

 
if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Run irradiance model using OOI profile input.')
#    parser.add_argument('--ooi_paramdict', help='A dictionary containing the necessary arguments to'\
#                                               ' download the ooi data. The keys to the dictionary'\
#                                               ' should be the following: [site, assembly, instrument,'\
#                                               ' method, start, stop]')
    parser.add_argument('--download', action='store_true', help='Download the data from OOI [This takes time], must specify file head')
    parser.add_argument('--load', action='store_true', help='Load the data from pickle, must specify file head.')
    parser.add_argument('--site_name', help='site name')
    parser.add_argument('--assembly', help= 'assembly')
    parser.add_argument('--method', help = 'method')
    parser.add_argument('--start', help='start')
    parser.add_argument('--stop', help='stop')
    

    
    parser.add_argument('--ooi_savefile_head', help='The name of the pickle file in which the ooi data is saved.'\
                                               'If this argument is given then the ooi data will be loaded'\
                                               ' from file and not downloaded from source.', 
                                               type=str, required=True)
    args = parser.parse_args()

    ## [The savefile names for each ooi data set.]
    spkir_savefile = args.ooi_savefile_head + '_spkir.p'
    flort_savefile = args.ooi_savefile_head + '_flort.p'
    optaa_savefile = args.ooi_savefile_head + '_optaa.p'

    ## [Download OOI data.]
    if args.download: 
        ## [Download the data.]
        ## [The Chla data.]
        flort_dat = Download_Data(args.site_name, args.assembly, 'FLORT', args.method, args.start, args.stop)
        ## [The spkir data.]
        spkir_dat = Download_Data(args.site_name, args.assembly, 'SPKIR', args.method, args.start, args.stop)
        ## [The OPTAA data.]
        optaa_dat = Download_Data(args.site_name, args.assembly, 'OPTAA', args.method, args.start, args.stop)
        ## [Save the ooi data into the pickle file.]
        pickle.dump(flort_dat, open(flort_savefile, 'wb'))
        pickle.dump(spkir_dat, open(spkir_savefile, 'wb'))
        pickle.dump(optaa_dat, open(optaa_savefile, 'wb'))

    ## [Load the data from pickle.]
    elif args.load:
        flort_dat = pickle.load(open(flort_savefile, 'rb'))
        spkir_dat = pickle.load(open(spkir_savefile, 'rb'))
        optaa_dat = pickle.load(open(optaa_savefile, 'rb'))

    ## [Get the chla profile lists.]
    depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs = ODF.Get_Flort_Profiles(flort_dat)

    ## [Get the spkir wavelengths.]
    spkir_wavelengths = np.array(ODF.Get_SPKIR_Wavelengths(spkir_dat.variables['spkir_abj_cspp_downwelling_vector']))

    ## [Run the irradiance model using the profiles, over all profiles.]
    N=100
    wavelengths = [410, 443, 486] #551, 638, 671]
    phy_types = ['HLPro', 'Cocco', 'Diat', 'Generic', 'Syn']
    PI = Param_Init()
    phy_type = phy_types[2]
    #irr_field, irr_field_ab = Run_Irradiance(N, wavelengths, phy_type, depth_profs, dt_profs, chla_profs, opt_bs_profs, optaa_dat)

    prof_index = 1
    ## [The index of the spkir wavelengths that most closely matches the given irradiance wavcelengths.]
    ## [This wavelength index should be automated soon.]
    #spkir_wavelength_index = [ODF.Get_Wavelength_Index(spkir_wavelengths, wavelengths[0])]
    #Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_wavelengths, spkir_wavelength_index, irr_field, irr_field_ab, spkir_dat, args.site_name, args.assembly, args.method)
    #Plot_Irr_OOI_Abs_Scat(PI, prof_index, wavelengths, N, phy_types, depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs, optaa_dat)

    #depth_ref = -10
    #Plot_OOI_Abs_Wavelength_Time(optaa_dat, depth_ref, args.start, args.stop, args.site_name, args.assembly, args.method)
    #Plot_OOI_Abs_Wavelength_Prof(optaa_dat, prof_index, args.start, args.stop, args.site_name, args.assembly, args.method)


    ## [Running the least square estimation of the ratio of phytoplankton.]
    A, b, x = Abs_Est_Species(PI, N, 443, prof_index, phy_types, depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs, optaa_dat)
