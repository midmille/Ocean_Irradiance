"""
est_ab_ooi.py
Author : Miles D. Miller University of California Santa Cruz 
Created : May 2, 2022

This file is for the estimation of the phytoplankton community that could compose the 
observed absorption and scattering in the water column.  

"""


## [Import statements.]
## [User Modules.]
import ooi_data_functions as ODF
from ocean_irradiance_module import Ocean_Irradiance as OI
import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module import Wavelength_To_RGB

## [External Modules.]
import numpy as np
import pickle 
from ooi_data_explorations.data_request import data_request
import matplotlib.pyplot as plt
import scipy


def Abs_Est_Species(PI, wavelengths, z_i, prof_index, phy_species, depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs, optaa_dat, ab='a', plot=False): 
    """
    This function is for the estimation of the ratio of different species that could compose the
    resulting observed absorption profile from OOI. The total observed phytoplankton concentration will be
    taken from the FLORT instrument and the observed absorption will be taken from the OPTAA instrument. 

    This estimation will be performed as a least squares problem where the absorption at different wavelengths will
    be considered the data points such that the system is overdetermined.

    This means the least squares will be implemented at a single location in z.

    The matrix problem is as follows: 

    Ax = y, where 

    ## [The effective absorption coefficient for each species at each wavelength]
    A = [[ a_phy1(lam_1), a_phy2(lam_1), a_phy3(lam_1), a_phy4(lam_1), ....], 
         [ a_phy1(lam_2), a_phy2(lam_2), a_phy3(lam_2), a_phy4(lam_2), ....],
         [    :  ,   :   ,  :    ,   :   , ....], 
         [    :  ,   :   ,  :    ,   :   , ....]]

    ## [The ratios of the different phytoplankton species, this is what we are solving for.]
    x = [ r_phy1, r_phy2, r_phy3, r_phy4, ....]^T 

    ## [The effective absorption given by observed OOI absorption divided by the 
        observed OOI chla and then subtracted by the dutkiewicz a_wat.]
    y = [a_ooi(lam_1), a_ooi(lam_2), a_ooi(lam_3), a_ooi(lam_4), a_ooi(lam_5) ... ]^T

    Note: The OOI absorption comes from the OPTAA instrument while the chla comes from the
    FLORT instrument. Thus it is necessary to interpolate them to the same grid. 

    Parameters
    ----------

    Returns 
    -------

    """

    ## [The number of species.]
    N_phy = len(phy_species)
    ## [The number of wavelengths.]
    N_lam = len(wavelengths)
    ## [The number of grid points for the averaged grid.]
    N_zavg = 100

    ## [The vertical point we are doing the wavelength based least square.]
    #z_i = 0

    ## [This step gets the OOI data.]
    ## [The optaa data.]
    optaa_depth_dat, optaa_dt_dat, optaa_abs_dat, optaa_scat, optaa_wavelength_dat = ODF.Get_Optaa_Dat(optaa_dat)

    dt_lbnd = dt_profs[prof_index].data[0] 
    dt_ubnd = dt_profs[prof_index].data[-1]
    prof_mask = ODF.OOI_Dat_Time_Prof_Mask(dt_lbnd, dt_ubnd, optaa_dt_dat)

    ## [This retreives the profiles]
    optaa_dt = optaa_dt_dat.data[prof_mask]
    z_ab = optaa_depth_dat.data[prof_mask]
    a_ooi = optaa_abs_dat.data[prof_mask]
    b_ooi = optaa_scat[prof_mask]
    optaa_wavelengths = optaa_wavelength_dat.data[prof_mask]

    ## [Get the arrays for the given prof_index.]
    ## [This is the chla grid and the grid to be used for everything.]
    z_chla = depth_profs[prof_index].data
    chla = chla_profs[prof_index].data
    cdom = cdom_profs[prof_index].data
    ## [Smoothing chla and cdom]
    z_chla_s, chla_s = ODF.Smooth_Profile_55(z_chla, chla)

    ## [The array of indexes that corresponds to desired wavelength.] 
    lam_is = np.zeros(N_lam, dtype=np.int8)
    ## [Loop over the wavelengths to get the corresponding ooi wavelength nearest neighbour.]
    for i in range(N_lam): 
        ## [Must get the wavelength index.]
        lam_is[i] = ODF.Get_Wavelength_Index(optaa_wavelengths[0,:], wavelengths[i])
        #lam_ooi = optaa_wavelengths[0,lam_i]

    ## [Construct the empty matrix system.]
    A = np.zeros((N_lam,N_phy))
    y = np.zeros(N_lam)

    ## [Loop over the wavelengths.]
    for i in range(N_lam): 

        ## [Absorption and scattering for current wavelength]
        a = np.squeeze(a_ooi[:,lam_is[i]])
        b = np.squeeze(b_ooi[:,lam_is[i]])
        
        ## [Smoothing the absorption and scattering for current wavelength]
        z_s, a_s = ODF.Smooth_Profile_55(z_ab, a)
        z_s, b_s = ODF.Smooth_Profile_55(z_ab, b)

        ## [Chla being interpolated to absorption grid.]
        chla_s_interp = np.interp(z_s, z_chla_s, chla_s)

        ## [Loop over the phytoplankton species.]
        for k in range(N_phy): 
            ## [Filling the A array, recall we are solving for a single z location. 
            ##  Each row corresponds to a wavelength and each column a different species.]
            ## [Absorption array.]
            if ab == 'a':
                A[i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[0] * chla_s_interp[z_i])
            if ab == 'b':
                A[i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[1] * chla_s_interp[z_i])

        ## [The rhs y, is simply the a_ooi minus the a_wat from Dut.]
        ## [The water absorption subtraction is commented out at the momment because it leads to negative.]

        ## [Get the cdom absorption]
        cdom_reflam = 412
        lam_cdom_i = ODF.Get_Wavelength_Index(optaa_wavelengths[0,:], cdom_reflam)
        CDOM = OI.CDOM_refa(z_s, 
                            a_s[z_i], 
                            optaa_wavelengths[0,lam_cdom_i],
                            optaa_wavelengths[0,lam_is[i]])
        a_cdom = CDOM.a

        ## WARNING 
        ## WARNING The abs should be removed and the negative absorption problem should be fixed, this is temp. 
        ## WARNING 
        ## [The cdom is subtracted.]
        if  ab == 'a':
            y[i] = abs(a_s[z_i]) - a_cdom # - abscat(wavelengths[i], 'water', C2chla='default')[0]
        if  ab == 'b':
            y[i] = abs(b_s[z_i]) 
        
    ## [Solving the least square system for the phy ratios.]
    #x = np.linalg.lstsq(A, b)
    x = scipy.optimize.nnls(A,y)[0]

    ## [Divide x by its two norm.]
    ## [x should be positive.]
#    x = x/ sum(x)

    #phy = OI.Phy(depth, chla_profs[k], esd(phy_type), 
    #             abscat(lam, phy_type, C2chla='default')[0], 
    #             abscat(lam, phy_type, C2chla='default')[1])

    #CDOM2C = 0
    #CDOM_dens = OI.CDOM_dens(depth, cdom_profs[k], CDOM2C, lam)

    if plot: 
               
        fig, ax = plt.subplots()

        y_fit = A@x

        ax.plot(wavelengths, y, 'o', label='OOI abs')
        ax.plot(wavelengths, y_fit, label='abs fit')
        ax.set_ylabel('Total Phy Abs [m^-1]')
        ax.set_xlabel('Wavelength [nm]')
        ax.grid()
        ax.legend()

        fig.show()

    return A, y, x 


def Plot_Abs_Est_LstSq_Sol(wavelengths, phy_species, A, y, x, ab='a'): 
    """
    This plots the resulting solution for the least squares absorption approximation. 
    
    The x-coordinate is the wavelength, and the y-coordinates is the effective absorption coefficient. 

    The absorption is plotted for each species used in the least squares estimation. The observed 
    total ooi effective absrroption is plotted and the resulting effective absorption coefficient
    from the approximated phytoplankton species ratio coeffeicients in plotted as well. 

    """
    
    ## [The number of phy species corresponds to the number of columns]
    N_phy = A.shape[1]

    y_est = A@x

    fig, ax = plt.subplots()

    if ab == 'a': 
        ab_str = 'Absorption'
    if ab == 'b': 
        ab_str = 'Scattering'

    ax.plot(wavelengths, y, 'k', label = f"OOI Total Effective {ab_str}")
    ax.plot(wavelengths, y_est, ':k', label = f"Least Square Approx {ab_str}")
    
    ## [Plot the species coefficients]
    for k in range(N_phy): 
        ax.plot(wavelengths, A[:, k], label = phy_species[k]) 

    ax.set_ylabel(f"Effective {ab_str}  [m^-1]")
    ax.set_xlabel("Wavelength [nm]")

    ax.grid()
    ax.legend()

    fig.show()

    return 


if __name__ == '__main__': 
    
    ## [The data request parameters.]
    site_name = "CE02SHSP"
    assembly = "profiler"
    method = "recovered_cspp"
    start = "2021-04-01"
    stop = "2021-04-30"

    ## [File name parameter.]
    ooi_savefile_head = "ooi_dat"

    ## [Data load flag, fale means load data from pickle.]
    download = False

    ## [Functional parameters.]
    N=100
    wavelengths = np.arange(425, 725, 25)
    #    phy_species = ['HLPro', 'LLPro', 'Cocco', 'Diat', 'Syn', 'Tricho', 'Lgeuk'] 
    phy_species = ['HLPro', 'Cocco', 'Diat', 'Syn'] 
    PI = Param_Init()
    cdom_reflam = 412.0
    prof_index = 0 
    z_i = -1
    ## [The absorption or scattering flag.]
    ab = 'b'


    ## [Download or load the data sets.]
    flort_dat, spkir_dat, optaa_dat = ODF.Download_OOI_Data(ooi_savefile_head, 
                                                            download, 
                                                            site_name, 
                                                            assembly, 
                                                            method, 
                                                            start, 
                                                            stop)
                                                            

    ## [Get the chla profile lists.]
    depth_profs, dt_profs, chla_profs, cdom_profs, opt_bs_profs = ODF.Get_Flort_Profiles(flort_dat)

    ## [Get the spkir wavelengths.]
    spkir_wavelengths = np.array(ODF.Get_SPKIR_Wavelengths(
                        spkir_dat.variables['spkir_abj_cspp_downwelling_vector']))


    ## [Running the least square estimation of the ratio of phytoplankton.]
    A, y, x = Abs_Est_Species(PI, 
                              wavelengths,
                                  z_i,
                                  prof_index,
                                  phy_species,
                                  depth_profs, 
                                  dt_profs, 
                                  chla_profs, 
                                  cdom_profs, 
                                  opt_bs_profs, 
                                  optaa_dat, 
                                  ab=ab,
                                  plot=False)


    Plot_Abs_Est_LstSq_Sol(wavelengths, phy_species, A, y, x, ab=ab)
