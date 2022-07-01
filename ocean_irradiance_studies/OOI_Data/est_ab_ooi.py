"""
est_ab_ooi.py
Author : Miles D. Miller University of California Santa Cruz 
Created : May 2, 2022

This file is for the estimation of the phytoplankton community that could compose the 
observed absorption and scattering in the water column.  

"""


## [Import statements.]
## [User Modules.]
import OOI_Data_Functions as ODF 
from ocean_irradiance_module import Ocean_Irradiance as OI
import ocean_irradiance_shubha.ocean_irradiance_shubha as OIS
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module import Wavelength_To_RGB
import plot_est_ab_ooi as PLOT

## [External Modules.]
import numpy as np
import pickle 
from ooi_data_explorations.data_request import data_request
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import cvxpy as cp


def Est_Spec_Lstsq(PI, wavelengths, depthz, phy_species, flort_prof, optaa_prof, ab='a', weight=None): 
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
    Nphy = len(phy_species)
    Nphyp1 = Nphy + 1
    Nphyp2 = Nphy + 2
    ## [The number of wavelengths.]
    Nlam = len(wavelengths)

    ## [Get the depth index for optaa.]
    zi_optaa = (abs(optaa_prof['depth'] - depthz)).argmin()
    ## [Get the depth index for flort.]
    zi_flort = (abs(flort_prof['depth'] - depthz)).argmin()

    ## [The reference wavelength for the cdom absorption.]
    cdom_reflam = 412

    ## [This retreives the profiles]
    optaa_dt = optaa_prof['time'].data
    optaa_z = optaa_prof['depth'].data
    wavelength_c = optaa_prof['wavelength_c'].data.transpose()
    wavelength_a = optaa_prof['wavelength_a'].data.transpose()
    optaa_c = optaa_prof['beam_attenuation'].data
    optaa_a = optaa_prof['optical_absorption'].data
    ## [interpolate the attenuation onto the absorption wavelength grid.]
#    optaa_interp = interpolate.interp1d(wavelength_c[0,:], optaa_, axis=1)
#    optaa_a = optaa_interp(wavelength_c[0,:])
    for k in range(len(optaa_z)): 
        optaa_c[k,:] = np.interp(wavelength_a[k,:], wavelength_c[k,:], optaa_c[k,:])

    print(optaa_a.shape)
    print(optaa_c.shape)
    ## [The scattering is simply the attenuation minus the absorption]
    optaa_b = optaa_c - optaa_a

    ## [The flort data.]
    flort_z = flort_prof['depth'].data
    flort_chla = flort_prof['fluorometric_chlorophyll_a'].data

    ## [The wieghting for the leats squares problem.]
    ## [If the weights are provided are argument, unpack and use.]
    if weight: 
        weighta, weightb = weight 
    ## [The default weight is the ratio of spectral mean absorption to spectral mean scattering.]
    else:
        ## [Since the weiht is a ratio of a/b, then only multiply scattering by weight.]
        weighta = 1 / np.nanmean(optaa_b[zi_optaa,:])
        weightb = 1 / np.nanmean(optaa_a[zi_optaa,:])

    ## [The array of indexes that corresponds to desired wavelength.] 
    lam_is = np.zeros(Nlam, dtype=np.int8)
    ## [Loop over the wavelengths to get the corresponding ooi wavelength nearest neighbour.]
    for i in range(Nlam): 
        ## [Must get the wavelength index.]
        lam_is[i] = ODF.Get_Wavelength_Index(wavelength_a[0,:], wavelengths[i])

    ## [The coupled scattering and absorption problem.]
    if ab == 'ab': 
        ## [Construct the empty matrix system.]
        ## [The Nphyp1 is to include the cdom estimation in the system.]
        A = np.zeros((2*Nlam,Nphyp2))
        y = np.zeros(2*Nlam)
    
        ## [Loop over the wavelengths.]
        for i in range(Nlam): 
    
            ## [Absorption and scattering for current wavelength]
            a = np.squeeze(optaa_a[:,lam_is[i]])
            b = np.squeeze(optaa_b[:,lam_is[i]])
            a_s = a
            b_s = b 
            z_s = optaa_z
         
            ## [Smoothing the absorption and scattering for current wavelength]
#            z_s, a_s = ODF.Smooth_Profile_55(optaa_z, a)
#            z_s, b_s = ODF.Smooth_Profile_55(optaa_z, b)
            
            ## [Chla being interpolated to absorption grid.]
            chla_s_interp = np.interp(optaa_z, flort_z, flort_chla )

            ## [Get the cdom absorption]
            lam_cdom_i = ODF.Get_Wavelength_Index(wavelength_a[0,:], cdom_reflam)
            CDOM = OI.CDOM_refa(z_s, 
                                a_s[zi_optaa], 
                                wavelength_a[0,lam_cdom_i],
                                wavelength_a[0,lam_is[i]])
            a_cdom = CDOM.a
    
            ## [Loop over the phytoplankton species.]
            for k in range(Nphy): 
                ## [Filling the A array, recall we are solving for a single z location. 
                ##  Each row corresponds to a wavelength and each column a different species.]
                ## [Absorption array.]
                ## [Multiply the absorption eqaution by weight a that is (1/(avg(absorptio(lam)))).]
                ## [Multiply the scattering equation by weight b. Which is (1/avg(scattering(lam))).] 
                A[2*i,k] = ((abscat(wavelengths[i], phy_species[k], C2chla='default')[0] * chla_s_interp[zi_optaa])) * weighta 

                A[2*i+1,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[1] * chla_s_interp[zi_optaa]) * weightb
            ## [Adding the cdom colummn.]
            A[2*i, Nphy] = a_cdom*weighta
            ## [The scattering for cdom is zero.]
            A[2*i+1, Nphy] = 0.0

            ## [Adding the detritus column.]
            A[2*i,Nphyp1] = abscat(wavelengths[i], 'detritus')[0] * weighta
            A[2*i+1,Nphyp1] = abscat(wavelengths[i], 'detritus')[1] * weighta
    
    
            ## WARNING The abs should be removed and the negative absorption problem should be fixed, this is temp. 
            ## [The cdom is subtracted.]
            y[2*i] = (abs(a_s[zi_optaa])) * weighta  
            y[2*i+1] = abs(b_s[zi_optaa]) * weightb


    ## [The uncoupled problem.]
    else:  
        ## [Construct the empty matrix system.]
        A = np.zeros((N_lam,N_phy))
        y = np.zeros(N_lam)
    
        ## [How m
    
        ## [Loop over the wavelengths.]
        for i in range(N_lam): 
    
            ## [Absorption and scattering for current wavelength]
            a = np.squeeze(optaa_a[:,lam_is[i]])
            b = np.squeeze(optaa_b[:,lam_is[i]])
            a_s = a
            b_s = b 
            z_s = optaa_z
 
            ## [Chla being interpolated to absorption grid.]
            chla_s_interp = np.interp(z_s, z_chla_s, chla_s)
    
            ## [Loop over the phytoplankton species.]
            for k in range(N_phy): 
                ## [Filling the A array, recall we are solving for a single z location. 
                ##  Each row corresponds to a wavelength and each column a different species.]
                ## [Absorption array.]
                if ab == 'a':
                    A[i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[0] * chla_s_interp[zi_optaa])
                if ab == 'b':
                    A[i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[1] * chla_s_interp[zi_optaa])
    
            ## [The rhs y, is simply the a_ooi minus the a_wat from Dut.]
            ## [The water absorption subtraction is commented out at the momment because it leads to negative.]
    
            ## [Get the cdom absorption]
            lam_cdom_i = ODF.Get_Wavelength_Index(wavelength_a[0,:], cdom_reflam)
            CDOM = OI.CDOM_refa(z_s, 
                                a_s[zi_optaa], 
                                wavelength_a[0,lam_cdom_i],
                                wavelength_a[0,lam_is[i]])
            a_cdom = CDOM.a
    
            ## WARNING 
            ## WARNING The abs should be removed and the negative absorption problem should be fixed, this is temp. 
            ## WARNING 
            ## [The cdom is subtracted.]
            if  ab == 'a':
                y[i] = abs(a_s[zi_optaa]) - a_cdom # - abscat(wavelengths[i], 'water', C2chla='default')[0]
            if  ab == 'b':
                y[i] = abs(b_s[zi_optaa]) 

    ## [Check for NaNs before solving.]
    ycheck = np.any(np.isnan(y))
    Acheck = np.any(np.isnan(A))
    if ycheck or Acheck: 
        y = y * np.nan
        A = A * np.nan
        x = np.zeros(N_phy) * np.nan
        print('NaN value found... Least squares estimation not implemented')
        return A, y, x

    ## [Solvving the problem as a convex optimization problem using the cvxpy.]
    ## [This is a least squares implementation fo the problem with constraints.]
    x = cp.Variable(Nphyp2)
    print('y', y)
    objective = cp.Minimize(cp.sum_squares(A@x - y))
    constraints = [0 <= x, cp.sum(x[:Nphy]) == 1.0]
    prob = cp.Problem(objective, constraints)
    result = prob.solve()
    x = x.value

    ## [If coupled then remove weight from returned array.]
    if ab == 'ab':
        A[0:2*Nlam:2, :] = A[0:2*Nlam:2, :] / weighta
        A[1:2*Nlam:2, :] = A[1:2*Nlam:2, :] / weightb
        y[0:2*Nlam:2] = y[0:2*Nlam:2] / weighta
        y[1:2*Nlam:2] = y[1:2*Nlam:2] / weightb
    ## [Solving the least square system for the phy ratios.]
    #x = np.linalg.lstsq(A, b)
#    x = scipy.optimize.nnls(A,y)[0]
#    x = scipy.optimize.lsq_linear(A,y, bounds=(0,1))

    ## [Divide x by its two norm.]
    ## [x should be positive.]
#    x = x/ sum(x)

    return A, y, x 


def Run_Lstsq_Time(N_profs, PI, wavelengths, depthz, phy_species, flort_profs, optaa_profs, ab='ab'): 
    """
    Run the least square problem over many profiles in time. The output from this function 
    can then be used to construct a hoffmuller diagram of the species ratios.
    """

    ## [The number of phytoplankton species.]
    N_phy = len(phy_species)

    ## [The exmpty x_profs array, which contains the lst sq species estimate for each profile
    ##  at the same depthz]
    x_profs = np.zeros((N_phy, N_profs))

    ## [The residuals of the fit for each.]
    residuals = np.zeros(N_profs)
    
    ## [Loop over the profiles.]
    for prof_index in range(N_profs):
        ## [Running the least square estimation of the ratio of phytoplankton.]
        A, y, x = Est_Spec_Lstsq(PI, 
                                wavelengths,
                                depthz,
                                phy_species,
                                flort_profs[prof_index], 
                                optaa_profs[prof_index],
                                ab=ab)
 
        ## [If the point is bad such that y is all zeros.]
        if sum(np.abs(y)) == 0.0: 
            x = np.NaN * x

        ## [Storing the result into x_profs array.]
        x_profs[:,prof_index] = x

        ## [The two norm of the difference between y_true and the apporx solution.]
        residuals[prof_index] = np.linalg.norm(A@x - y, ord=2) / np.linalg.norm(y, ord=2)

    return x_profs, residuals


def Run_Lstsq_Depth(N_depths, PI, wavelengths, phy_species, flort_prof, optaa_prof, ab='ab'): 
    """
    Run the least square problem over a sinlge profile for all optaa depths. The output from this function 
    can then be used to construct a hoffmuller diagram of the species ratios.
    """

    ## [The number of phytoplankton species.]
    N_phy = len(phy_species)

    ## [Get the depths from the optaa data set.]
    depths = optaa_prof['depth'].data

    ## [The exmpty x_profs array, which contains the lst sq species estimate for each profile
    ##  at the same depthz]
    x_profs = np.zeros((N_phy, N_depths))

    ## [The residuals of the fit for each.]
    residuals = np.zeros(N_depths)
    
    ## [Loop over the profiles.]
    for k in range(N_depths):
        ## [Running the least square estimation of the ratio of phytoplankton.]
        A, y, x = Est_Spec_Lstsq(PI, 
                                wavelengths,
                                depths[k],
                                phy_species,
                                flort_prof, 
                                optaa_prof,
                                ab=ab)
 
        ## [If the point is bad such that y is all zeros.]
        if sum(np.abs(y)) == 0.0: 
            x = np.NaN * x

        ## [Storing the result into x_profs array.]
        x_profs[:,k] = x

        ## [The two norm of the difference between y_true and the apporx solution.]
        residuals[k] = np.linalg.norm(A@x - y, ord=2) / np.linalg.norm(y, ord=2)

    return x_profs, residuals


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

    ## [Functional parameters.]
    N=100
    wavelengths = np.arange(425, 725, 25)
#    phy_species = ['HLPro', 'LLPro', 'Cocco', 'Diat', 'Syn', 'Lgeuk'] 
    phy_species = [ 'Cocco', 'Diat', 'Syn', 'Lgeuk'] 
#    phy_species = ['HLPro', 'Cocco', 'Diat', 'Syn'] 
    PI = Param_Init()
    cdom_reflam = 412.0
    prof_index = 1
    depthz = -9
    ## [The absorption or scattering flag.]
    ab = 'ab'
    ## [The number of profiles for the Hoffmuller.]
    N_profs = 50
    ## [Number of depths for the depth Hoffmuller.]
    N_depths = 270
#    N_depths = 20

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

    ## [Get the spkir wavelengths.]
    spkir_wavelengths = np.array(ODF.Get_SPKIR_Wavelengths(
                        spkir_dat.variables['spkir_abj_cspp_downwelling_vector']))

    ## [Running the least square estimation of the ratio of phytoplankton.]
    A, y, x = Est_Spec_Lstsq(PI, 
                             wavelengths,
                             depthz,
                            phy_species,
                             flort_prof, 
                             optaa_prof,
                             ab=ab)
    PLOT.Plot_Abs_Est_LstSq_Sol(wavelengths, optaa_prof, phy_species, A, y, x, color_dict, ab=ab)

    ## [Running the least squares problem over many profiles.]
#    x_profs, residuals = Run_Lstsq_Time(N_profs, 
#                                       PI, 
#                                        wavelengths, 
#                                        depthz, 
#                                        phy_species, 
#                                        flort_profs, 
#                                        optaa_profs,
#                                        ab=ab)
#    PLOT.Plot_Spec_Lstsq_Hoffmuller_Time(N_profs, depthz, flort_profs, optaa_profs, phy_species, x_profs, color_dict, plot_residuals=True, residuals=residuals)

    ## [Running the least squares problem over many depths in single profile.]
#    x_profs, residuals = Run_Lstsq_Depth(N_depths, 
#                                        PI, 
#                                        wavelengths, 
#                                        phy_species, 
#                                      flort_prof, 
#                                        optaa_prof,
#                                       ab=ab)
#    PLOT.Plot_Spec_Lstsq_Hoffmuller_Depth(N_depths, flort_prof, optaa_prof, phy_species, x_profs, color_dict, plot_residuals=True, residuals=residuals)
