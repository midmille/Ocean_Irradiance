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
import OOI_Data_Functions as ODF 
from ocean_irradiance_module import Ocean_Irradiance as OI
from ocean_irradiance_module import Ocean_Irradiance_ROMS as OIR
import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module import Wavelength_To_RGB
import ocean_irradiance_visualization.Plot_Comparison as PC
import plymouth_data

## [External Modules]
import numpy as np
import pickle
from ooi_data_explorations.data_request import data_request
import matplotlib.pyplot as plt
import scipy
import geopy.distance


def OOI_Abs_Scat(optaa_prof, lam, smooth=True):
    """
    This function gets the OOI absorption and scattering for a given profile index and data sets.
    """

    optaa_z = optaa_prof['depth'].data
    wavelength_c = optaa_prof['wavelength_c'].data.transpose()
    wavelength_a = optaa_prof['wavelength_a'].data.transpose()
    optaa_c = optaa_prof['beam_attenuation'].data
    optaa_a = optaa_prof['optical_absorption'].data
    for k in range(len(optaa_z)): 
        optaa_c[k,:] = np.interp(wavelength_a[k,:], wavelength_c[k,:], optaa_c[k,:])

    ## [The scattering is simply the attenuation minus the absorption]
    optaa_b = optaa_c - optaa_a

    ## [Must get the wavelength index.]
    lam_i = ODF.Get_Wavelength_Index(wavelength_a[0,:], lam)
    lam_ooi = wavelength_a[0,lam_i]

    ##[Labeling the ooi z-grids and a, b_b.]
    z = optaa_z
    ## [The squeese makes the array 1-D.]
    ## [Removing the water.]
    a = np.squeeze(optaa_a[:, lam_i])
    b = np.squeeze(optaa_b[:, lam_i])

    ## [Smoothing the data using the 55 Smoothing algorithm.] 
#    z_s, a_s = ODF.Smooth_Profile_55(z, a)
#    z_s, b_s = ODF.Smooth_Profile_55(z, b)
    z_s = z
    a_s = a
    b_s = b

    return z_s, a_s, b_s, lam_ooi


def Irr_OOI_Abs_Scat(PI, N, lam, phy_type, flort_prof, optaa_prof,  cdom_reflam): 
    """
    This function compares the absorption and backscatter as calculated by the irradiance model with Dutkiewicz
    coefficients to the values given by OO.
    """
    
    chla = flort_prof['fluorometric_chlorophyll_a'].data
    flort_z = flort_prof['depth'].data
    phy = OI.Phy(flort_z, chla, esd(phy_type), 
                 abscat(lam, phy_type, C2chla='default')[0], 
                 abscat(lam, phy_type, C2chla='default')[1])

    ## [Get the z_a, a, and lam_ooi for the CDOM_refa object.]
    z_a, cdom_refa, b, lam_ooi = OOI_Abs_Scat(optaa_prof, cdom_reflam, smooth=True)

    ## [Assumes that all absorption at the smallest wavelength is due to CDOM.]
    CDOM = OI.CDOM_refa(z_a, cdom_refa, cdom_reflam, lam, fraca=1.0)

    print('CDOM', CDOM)
            
    z_irr, a_irr, b_irr, b_b_irr  = OIS.ocean_irradiance_two_stream_ab(flort_z[0], 
                                                                   abscat(lam, 'water'), 
                                                                   N,
                                                                   phy=phy, 
                                                                   CDOM_refa=CDOM, 
                                                                   pt1_perc_zbot=False, 
                                                                   pt1_perc_phy=False)

    ## [Get the ooi absorption and scattering, along with their corresponding grids.]
    z_ooi, a_ooi, b_ooi, lam_ooi = OOI_Abs_Scat(optaa_prof, lam, smooth=True)


    return  z_irr, a_irr, b_irr, z_ooi, a_ooi, b_ooi, lam_ooi


def Run_Irradiance(PI, N, wavelengths, spkir_wavelengths, phy_type, flort_profs, optaa_profs, spkir_profs, cdom_reflam):
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
#    N_profs = len(flort_profs)
    N_profs = 30

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
            flort_prof = flort_profs[k]
            optaa_prof = optaa_profs[k]

            chla = flort_prof['fluorometric_chlorophyll_a'].data
            flort_z = flort_prof['depth'].data
            phy = OI.Phy(flort_z, chla, esd(phy_type), 
                        abscat(lam, phy_type, C2chla='default')[0], 
                        abscat(lam, phy_type, C2chla='default')[1])
    
            ## [Get the z_a, a, and lam_ooi for the CDOM_refa object.]
            z_a, cdom_refa, b, lam_ooi = OOI_Abs_Scat(optaa_prof, cdom_reflam, smooth=True)

            ## [Assumes that all absorption at the smallest wavelength is due to CDOM.]
            CDOM = OI.CDOM_refa(z_a, cdom_refa, cdom_reflam, lam, fraca=1.0)

            ## [The detritus parameter.]
            det_val = 0
            det = OI.Det(flort_z, np.full(len(flort_z), det_val), 
                         abscat(lam, 'detritus')[0], 
                         abscat(lam, 'detritus')[1], 
                         abscat(lam, 'detritus')[2])
            
            z_irr, a_irr, b_irr, b_b_irr  = OIS.ocean_irradiance_two_stream_ab(flort_z[0], 
                                                                            abscat(lam, 'water'), 
                                                                            N,
                                                                            phy=phy, 
                                                                            CDOM_refa=CDOM, 
                                                                            det = det,
                                                                            pt1_perc_zbot=False, 
                                                                            pt1_perc_phy=False)

            ## [This solves for the irradiance solution using Dutkiewicz coefficients and 
            ##  CDOM.]
            Ed0 = spkir_profs[k]['spkir_abj_cspp_downwelling_vector'][-1, ODF.Get_Wavelength_Index(spkir_wavelengths, lam)]
            ocean_irr_sol = OIS.ocean_irradiance_shubha_ooi(flort_z[0], 
                                                        Ed0, 
                                                        PI.coefficients, 
                                                        z_irr, 
                                                        a_irr, 
                                                        z_irr, 
                                                        b_b_irr,
                                                        N=N, 
                                                        pt1_perc_zbot=False)

            ## [Storing output into array.]
            irr_arr[:,k,0] = ocean_irr_sol[0]
            irr_arr[:,k,2] = ocean_irr_sol[1]
            irr_arr[:,k,3] = ocean_irr_sol[2] 

#            ## [Now using the OOI absorption and scattering.]
#            ## [The optaa data.]
#            optaa_depth_dat = optaa_dat.variables['depth'].data
#            optaa_dt_dat = optaa_dat.variables['time'].data
#            optaa_abs_dat = optaa_dat.variables['optical_absorption'].data
#            optaa_wavelength_dat = optaa_dat.variables['wavelength_a'].data
#
#            dt_lbnd = dt_profs[k].data[0] 
#            dt_ubnd = dt_profs[k].data[-1]
#            prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, optaa_dt_dat)
#
#            optaa_dt_prof = optaa_dt_dat[prof_mask]
#            optaa_depth_prof = optaa_depth_dat[prof_mask]
#            optaa_abs_prof = optaa_abs_dat[prof_mask]
#            optaa_wavelengths = optaa_wavelength_dat[prof_mask]
#
#            optaa_depth = -np.squeeze(optaa_depth_prof)
#            ## [Must get the wavelength index.]
#            lam_i = ODF.Get_Wavelength_Index(optaa_wavelengths[0,:], lam)
#
#            ## [The squeese makes the array 1-D.]
#            a = np.squeeze(optaa_abs_prof[:, lam_i])
#            b = opt_bs_profs[k]
            
            ## [Get the absorbtion and scattering for OOI.]
            z_ooi, a_ooi, b_ooi, lam_ooi = OOI_Abs_Scat(optaa_prof, lam) 

            ## [Change the scattering into backscattering.]
            esd_v = phy.esd
            ## [backscatter ratio.]
            bb_r = OI.Backscatter_Ratio(esd_v)

            ## [Add water into the absorption/scattering.]
            a_ooi = a_ooi + abscat(lam, 'water')[0]
            bb_wat = 0.551 * abscat(lam, 'water')[1]
            b_ooi = b_ooi*bb_r + bb_wat

            ## [Now running the irr model using abs and scat from ooi.]
            ocean_irr_sol_ab = OIS.ocean_irradiance_shubha_ooi(z_ooi[0], 
                                                              Ed0, 
                                                              PI.coefficients, 
                                                              z_ooi,  
                                                              a_ooi, 
                                                              z_ooi, 
                                                              b_ooi,
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


def Get_CCI_Data(pml_url, flort_dat): 
    """
    This gets the CCI data from the PML servers that closest corresponds to the OOI data time and location.
    """

    ## [The year limits from the flort data set.]
    year_lims = [ int(str(flort_dat['time'].data.astype('datetime64[Y]')[0])), 
                  int(str(flort_dat['time'].data.astype('datetime64[Y]')[-1]))]

    ## [The julian date limits.]
    jd_lims = [str(flort_dat['time'].data.astype('datetime64[D]')[0]), 
               str(flort_dat['time'].data.astype('datetime64[D]')[-1])]
    for k, d in enumerate(jd_lims): 
        jd_lims[k] = plymouth_data.Date_to_Julian_Date(d)

    ## [The lattitude limits.]
    lat = flort_dat.attrs['geospatial_lat_min']
    lat_lims = [lat[0] + 0.3, lat[0] - 0.3]

    ## [The longitude limits.]
    lon = flort_dat.attrs['geospatial_lon_min']
    lon_lims = [lon[0] - 0.5, lon[0] + 0.5]

    ## [Download the pml cci data set.]
    cci_dat = plymouth_data.Download_PML_Data(pml_url, year_lims, jd_lims, lat_lims, lon_lims)

    return cci_dat


def Get_OOI_CCI_Match(cci_ds, flort_dat, flort_profs): 
    """
    This function is for finding the CCI nearest neighbour to the OOI mooring.
    """

    ## [The latitude and longitude coordinates for the OOI platform.]
    ooi_lat = flort_dat.attrs['geospatial_lat_min'][0]
    ooi_lon = flort_dat.attrs['geospatial_lon_min'][0]

    ## [The latitured and longitude of cci data.]
    cci_lat = cci_ds.variables['lat'][:]
    cci_lon = cci_ds.variables['lon'][:]

    ## [The cci time array.]
    time = cci_ds.variables['time'][:]

    ## [The CCI chla mask.]
    chla = cci_ds.variables['chlor_a'][:]
    mask = chla.mask

    ## [Init the return array.]
    ## [The four is for the resulting indexes i.e.
    ##  res[:,0] = dti, res[:,1] = loni, res[:,2] = lati, res[:,3] = distnn
    res = np.zeros((len(flort_profs), 4)) 

    ## [loop over the OOI profiles.]  
    for k, prof in enumerate(flort_profs): 

        ## [First find the CCI data fro the given day.]
        ## [Convert the day of the flort profile to julian day.]
        ooi_dt = str(prof['time'].data.astype('datetime64[D]')[0])
        ## [This concatinates the year with julain date.]
        ooi_jd = plymouth_data.Date_to_Julian_Date(ooi_dt)
        ## [Conver the CCI data from days since 1970 to julian day, year.]
        cci_jd, cci_year = plymouth_data.Days_To_Julian_Date('01/01/1970',time)
        ## [The date time index.]
        dti = np.where(cci_jd == ooi_jd)[0]
        
        ## [Loop over the lon.]
        cnt = 0
        ## [Init values s.t. if the mask is false or chla>1e30 then NaN.]
        jnn = 0
        inn = 0
        distnn = np.NaN

        ## [print dti]
        print(dti)
        for i, lon in enumerate(cci_lon): 
            ## [Loop over the lat.]
            for j, lat in enumerate(cci_lat): 
                ## [check the existence of chla.]
                if (mask[dti, j, i] == True) and (chla[dti, j, i] < 1e30): 
                    
                    ## [Calculate the distance from the OOi morring.]
                    dist = geopy.distance.distance((ooi_lat, ooi_lon), (lat, lon)).kilometers

                    ## [Init for the first loop.]
                    if cnt == 0: 
                        ## [Nearest neighbour distance.]
                        distnn = dist
                        inn =  i 
                        jnn = j
                    elif dist < distnn: 
                        ## [Nearest neighbour distance.]
                        distnn = dist
                        inn =  i 
                        jnn = j

                    cnt += 1
 
        ## [Write the result.]
        res[k,:] = np.array([dti, jnn, inn, distnn])


    return res


def Calc_Chla_Irr(prof_index, irr_field): 
    """
    This function calculates the chla at the surface using the OCx algorithm and RRs ratios. It takes the 
    irradiance field as an argument and returns the chla.

    ASSUMES TWO STREAM IRRADIANCE
    """

    ## [The blue wavelength.]
    lamb = 443
    ## [The green wavelength.]
    lamg = 551

    rrsb = OIR.R_RS(irr_field[lamb][-1,prof_index,0], 0, irr_field[lamb][-1,prof_index,2])
    rrsg = OIR.R_RS(irr_field[lamg][-1,prof_index,0], 0, irr_field[lamg][-1,prof_index,2])

    chla = OIR.OCx_alg(rrsb, rrsg)


    return chla


def Comp_OOI_CCI_Irr(PI, N, wavelengths, spkir_wavelengths, phy_species, cci_ds, flort_dat, flort_profs, optaa_profs, spkir_profs, cdom_reflam, plot=True): 
    """
    This function compares the ooi flourometric chla, the irradiance chla with dutkiewicz absorption/scattering,
    the irradiance chla with ooi abs/scat, and the chla from cci satellite.
    """
    
    ## [The numer of OOI profiles.]
#    Nprofs = len(flort_profs)
    Nprofs = 30
    ## [The number of phy_species.]
    Nphy = len(phy_species)

    ## [The empty chla arrays.]
    ## [This one is the Dut. ab coefficients.]
    ooi_chla = np.zeros((Nprofs, Nphy))
    ## [This one is ab from ooi.]
    ooi_chla_ab = np.zeros((Nprofs, Nphy))
    ## [Fluorometric chla.]
    flort_chla = np.zeros((Nprofs, Nphy))
    ## [The CCI chlor_A]
    cci_chla = np.zeros((Nprofs, Nphy))


    ## [Loop the phy_species.]
    for ip, phy_type in enumerate(phy_species):     ## [Calculate the irradiance feilds.]

        irr_field, irr_field_ab = Run_Irradiance(PI, N, wavelengths, spkir_wavelengths, phy_type, flort_profs, optaa_profs, spkir_profs, cdom_reflam)

        ## [Get the nearest neighbour CCI indexes.]
        res = Get_OOI_CCI_Match(cci_ds, flort_dat, flort_profs)
    
        ## [Loop over the time coordinates of the cci.]
        for k in range(Nprofs): 
            ## [The indexes at the given oor profile.] 
            dti = int(res[k,0])
            jnn = int(res[k,1])
            inn = int(res[k,2])
            dist =  res[k,3]
    
            ## [If the CCI data is bad for the entire search square ignore it.]
            if dist == np.NaN: 
                ooi_chla[k,ip] = np.NaN
                ooi_chla_ab[k,ip] = np.NaN
                flort_chla[k,ip] = np.NaN
                cci_chla[k,ip] = np.NaN
    
            else:
                ## [Calculate the chla from irradiance.]
                ooi_chla[k,ip] = Calc_Chla_Irr(k, irr_field)
                ooi_chla_ab[k,ip] = Calc_Chla_Irr(k, irr_field_ab)
    
                ## [Get the fluorometric chla from data.]
                flort_chla[k,ip] = flort_profs[k].variables['fluorometric_chlorophyll_a'].data[-1]
    
                ## [Get the chla from the cci_ds.]
                cci_chla[k,ip] = cci_ds.variables['chlor_a'][:][dti, jnn, inn]

                if ooi_chla[k,ip] > 100: 
                    print('ATTENTION BAD profile', k)
                    irr_field_save = irr_field
                   


    if plot:  
        
        fig, axs = plt.subplots(ncols=3, nrows=2)

        for ip, phy_type in enumerate(phy_species):

            title = 'CCI Chla against OOI \n Irradiance Chla Dut. ab Coefficients'
            label = phy_type
            xlabel = 'CCI Chla' 
            ylabel = 'OOI Chla Irr'
            PC.Plot_Comparison(axs[0,0], cci_chla[:,ip], ooi_chla[:,ip], title, label, xlabel, ylabel)
    
            title = 'CCI Chla against OOI \n Irradiance Derived Chla OOI ab Coefficients'
            xlabel = 'CCI Chla' 
            ylabel = 'OOI Chla Irr'
            PC.Plot_Comparison(axs[0,1], cci_chla[:,ip], ooi_chla_ab[:,ip], title, label, xlabel, ylabel)
            
            title = 'CCI Chla against OOI Fluorometric Chla'
            xlabel = 'CCI Chla' 
            ylabel = 'OOI Chla Flort'
            PC.Plot_Comparison(axs[0,2], cci_chla[:,ip], flort_chla[:,ip], title, label, xlabel, ylabel)
     
            title = 'OOI fluorometric chla against OOI\n  Irradiance Derived Chla Dut. ab Coefficients'
            xlabel = 'OOI Chla Flort' 
            ylabel = 'OOI Chla Irr'
            PC.Plot_Comparison(axs[1,0], flort_chla[:,ip], ooi_chla[:,ip], title, label, xlabel, ylabel)
    
            title = 'OOI fluorometric chla against OOI \n Irradiance Derived Chla OOI ab Coefficients'
            xlabel = 'OOI Chla Flort' 
            ylabel = 'OOI Chla Irr ab'
            PC.Plot_Comparison(axs[1,1], flort_chla[:,ip], ooi_chla_ab[:,ip], title, label, xlabel, ylabel)

            title = 'OOI Irradiance Derived Chla OOI ab Coefficients'
            xlabel = 'OOI Chla Irr ab' 
            ylabel = 'OOI Chla Irr'
            PC.Plot_Comparison(axs[1,2], ooi_chla_ab[:,ip], ooi_chla[:,ip], title, label, xlabel, ylabel)

        for ax in axs.flatten(): 
            ax.legend()

        fig.show()

        
    return ooi_chla, ooi_chla_ab, cci_chla



def Plot_Irr_OOI_Abs_Scat(PI, wavelengths, N, phy_types, flort_prof, optaa_prof, cdom_reflam): 
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
            z_irr, a_irr, b_irr, z_ooi, a_ooi, b_ooi, lam_ooi = Irr_OOI_Abs_Scat(PI, N, lam, phy_type, flort_prof, optaa_prof, cdom_reflam)

            ## [Adding water absorption to OOI absorption.]
#            a_ooi = a_ooi + abscat(lam, 'water')[0]
#            b_ooi = b_ooi + abscat(lam, 'water')[1]
            ## [Removing water from irr.]
            a_irr = a_irr - abscat(lam, 'water')[0]
            b_irr = b_irr - abscat(lam, 'water')[1]


            ## [Plotting Absorption comparison.]
            ## [Only plot the ooi abs, scat for i==0, since they are indp. of wavelength.]
            if i == 0: 
                axa.plot(a_ooi, z_ooi, ':', label=f'Abs OOI OPTAA {lam_ooi}')
            axa.plot(a_irr, z_irr, label=f'Abs Irr {phy_type}')
            axa.set_title(f'Absorption {lam}')
            axa.set_ylabel('Z[m]')
            axa.set_xlabel('Abs [m^-1]')
            axa.legend()
            axa.grid()
    
            ## [Plotting the scattering comp.]
            ## [Only plot the ooi abs, scat for i==0, since they are indp. of wavelength.]
            if i == 0: 
                axb.plot(b_ooi, z_ooi, ':', label=f'Scat OOI FLORT {lam_ooi}')
            axb.plot(b_irr, z_irr, label=f'Scat Irr {phy_type}')
            axb.set_title(f'Total Scattering {lam}')
            axb.set_ylabel('Z[m]')
            axb.set_xlabel('Scat [m^-1]')
            axb.legend()
            axb.grid()
 
    
    figa.show()
    figb.show()

    return 


def Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_prof, spkir_wavelengths, irr_field, irr_field_ab, site, method): 

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
    spkir_depth = - spkir_dat['depth'].data
    spkir_dt = spkir_dat['time'].data
    spkir = spkir_dat['spkir_abj_cspp_downwelling_vector'].data

    fig, ax = plt.subplots()
    ## [Get the colors such that they match the wavelengths.] 
    colors = [Wavelength_To_RGB.wavelength_to_rgb(wavelength) for wavelength in wavelengths]
    
    ## [Loop over the irradiance wavelengths.]
    for k, lam in enumerate(wavelengths): 
        ## [Get the spkir wavelength index for lam.]
        lam_i = ODF.Get_Wavelength_Index(spkir_wavelengths, lam)
        lam_spkir = spkir_wavelengths[lam_i]

        ## [Make the 1m avg grid.]
        depth_avg = np.arange(spkir_depth[0], spkir_depth[-1], 1)
        spkir_avg = ODF.Grid_Average_Profile(spkir_depth, spkir[:,lam_i], depth_avg)
        
        ## [Surface Ed0.]
        print('surface Ed0+Es0', spkir_avg[-1])

        ## [irr arrays for lam.]
        irr_arr = irr_field[lam]
        irr_arr_ab = irr_field_ab[lam]
        ## [ Plotting the downward direct profile. 0 is downward irradiance, 3 is depth.]
        ## [Also must multiply by surface value of spkir prof, since irr_arr is normalized.]
        ax.plot(spkir_avg[-1] * irr_arr[:, prof_index, 0], irr_arr[:, prof_index, 3], ':', label=f'Irr {lam}', color=colors[k])
        ax.plot(spkir_avg[-1] * irr_arr_ab[:, prof_index, 0], irr_arr_ab[:, prof_index, 3], '-', label=f'Irr ooi ab {lam}', color=colors[k])

        ## [Plotting the spkir profile.]
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
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: SPKIR')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.legend(title='Wavelengths [nm]')
    ax.grid()

    fig.show()

    return 

def Plot_OOI_Abs_Wavelength_Time(optaa_dat, depth_ref, start, stop, site, assembly, method):
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
    ax.text(txt_x, txt_y, f'SITE: {site}')   
    ax.text(txt_x, txt_y+txt_dz, f'ASSEMBLY: {assembly}')   
    ax.text(txt_x, txt_y+2*txt_dz, f'INSTRUMENT: OPTAA')   
    ax.text(txt_x, txt_y+3*txt_dz, f'METHOD: {method}')   

    ax.grid()
    
    fig.show()

    return
    

def Plot_OOI_Abs_Wavelength_Prof(optaa_dat, prof_index, start, stop, site, assembly, method):
    """
    This plot creates a plot similiar to the one Chris Wingard sent in his email on 4/9/2022. 
    The main difference between this function and Plot_OOI_Abs_Wavelength_Time() is that this function
    plots the absorption against wavelength for all depths for a single profile. 
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

    return


def Plot_OOI_CCI_Loc(cci_ds, flort_dat, flort_profs): 
    """
    This function plots the location comparisons of the CCI and OOI data.
    """

    ## [The latitude and longitude coordinates for the OOI platform.]
    ooi_lat = flort_dat.attrs['geospatial_lat_min'][0]
    ooi_lon = flort_dat.attrs['geospatial_lon_min'][0]

    ## [The latitured and longitude of cci data.]
    cci_lat = cci_ds.variables['lat'][:]
    cci_lon = cci_ds.variables['lon'][:]

    ## [Plotting Result]
    fig, ax = plt.subplots()

    ## [plot the location of the OOI mooring.]
    ax.plot(ooi_lon, ooi_lat, 'k*')

    ## [Finding the cci indexing of the nearest neighbour satellite point to OI mooring.]
    res = Get_OOI_CCI_Match(cci_ds, flort_dat, flort_profs)

    ## [Loop the flort profiles.]
    for k, prof in enumerate(flort_profs): 
        dti = int(res[k,0])
        jnn = int(res[k,1])
        inn = int(res[k,2])
        dist =  res[k,3]
        print("Dist:", dist)
        ## [Plot the OOI nearest neighbour.]
        ax.plot(cci_lon[inn], cci_lat[jnn], 'bo')

    fig.show()

        
    return 


if __name__ == '__main__':

    ## [The data request parameters.]
    ## [Data load flag, false means load data from pickle.]
    download = False
    savefile_head = "ooi_data/ooi_dat"
    site = "CE02SHSP"
    node = 'SP001'
    method = 'recovered_cspp'
    deployment = 15
    ## [Makes the OOI data into profiles.]
    profiled = True
    ## [Processes the optaa data according to Chris Winegards script.]
    optaa_process_raw = True

    ## [CCI sattelite data parameters.]
    pml_url = 'https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v5.0-DAILY?lat[0:1:0],lon[0:1:0],time[0:1:0],Rrs_443[0:1:0][0:1:0][0:1:0],Rrs_560[0:1:0][0:1:0][0:1:0],chlor_a[0:1:0][0:1:0][0:1:0]' 

    ## [Functional parameters.]
    ## [The number of levels for irradiance run.]
    N=1000
#    wavelengths = np.arange(425, 725, 25)
    wavelengths = [443, 551]
#    phy_species = ['HLPro', 'LLPro', 'Cocco', 'Diat', 'Syn', 'Lgeuk'] 
    phy_species = [ 'Syn'] 
    PI = Param_Init()
    cdom_reflam = 412.0
    prof_index = 0
    depthz = 5

    ## [Get color list for each species of phy.]
    color_dict = ODF.Get_Phy_Cmap_Dict()

    ## [Download or load the data sets.]
    ooi_data = ODF.Download_OOI_Data(savefile_head, 
                                     download, 
                                     site,
                                     node,
                                     method, 
                                     deployment,
                                     profiled = profiled, 
                                     optaa_process_raw = optaa_process_raw)

    flort_dat, spkir_dat, optaa_dat, flort_profs, spkir_profs, optaa_profs = ooi_data

    ## [Index the profiles for the given index.]
    flort_prof = flort_profs[prof_index]
    optaa_prof = optaa_profs[prof_index]
    spkir_prof = spkir_profs[prof_index]

    ## [Download the relevant CCI data set.]
    cci_ds = Get_CCI_Data(pml_url, flort_dat)

    ## [Get the spkir wavelengths.]
    spkir_wavelengths = np.array(ODF.Get_SPKIR_Wavelengths(
                        spkir_dat.variables['spkir_abj_cspp_downwelling_vector']))

    ## [Run the irradiance model using the profiles, over all profiles.]
    phy_type = 'Syn'
    #irr_field, irr_field_ab = Run_Irradiance(PI, N, wavelengths, phy_type, flort_profs, optaa_profs, cdom_reflam)
    ## [Plot the resulting irradiance profiles.]
#    Plot_Irraddiance_SPKIR(prof_index, wavelengths, spkir_prof, spkir_wavelengths, irr_field, irr_field_ab, site, method)

#    Plot_Irr_OOI_Abs_Scat(PI, wavelengths, N, phy_species, flort_prof, optaa_prof, cdom_reflam)

    ooi_chla, ooi_chla_ab, cci_chla = Comp_OOI_CCI_Irr(PI, N, wavelengths, spkir_wavelengths, phy_species, cci_ds, flort_dat, flort_profs, optaa_profs, spkir_profs, cdom_reflam)
