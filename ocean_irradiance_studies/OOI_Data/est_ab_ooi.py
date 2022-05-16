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

## [External Modules.]
import numpy as np
import pickle 
from ooi_data_explorations.data_request import data_request
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import cvxpy as cp


def Est_Spec_Lstsq(PI, wavelengths, z_i, phy_species, flort_prof, optaa_prof, ab='a', plot=False): 
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

    ## [The reference wavelength for the cdom absorption.]
    cdom_reflam = 412

    ## [This retreives the profiles]
    optaa_dt = optaa_prof['time'].data
    optaa_z = optaa_prof['depth'].data
    wavelength_c = optaa_prof['wavelength_c'].data
    wavelength_a = optaa_prof['wavelength_a'].data
    optaa_c = optaa_prof['beam_attenuation'].data
    optaa_a = optaa_prof['optical_absorption'].data
    ## [interpolate the attenuation onto the absorption wavelength grid.]
    optaa_interp = interpolate.interp1d(wavelength_a[0,:], optaa_a, axis=1)
    optaa_c = optaa_interp(wavelength_c[0,:])

    print(optaa_a.shape)
    print(optaa_c.shape)
    ## [The scattering is simply the attenuation minus the absorption]
    optaa_b = optaa_c - optaa_a

    ## [The flort data.]
    flort_z = flort_prof['depth'].data
    flort_chla = flort_prof['fluorometric_chlorophyll_a'].data
    z_chla_s, chla_s = ODF.Smooth_Profile_55(flort_z, flort_chla)


    ## [The array of indexes that corresponds to desired wavelength.] 
    lam_is = np.zeros(N_lam, dtype=np.int8)
    ## [Loop over the wavelengths to get the corresponding ooi wavelength nearest neighbour.]
    for i in range(N_lam): 
        ## [Must get the wavelength index.]
        lam_is[i] = ODF.Get_Wavelength_Index(wavelength_a[0,:], wavelengths[i])

    ## [The coupled scattering and absorption problem.]
    if ab == 'ab': 
        ## [Construct the empty matrix system.]
        A = np.zeros((2*N_lam,N_phy))
        y = np.zeros(2*N_lam)
    
        ## [Loop over the wavelengths.]
        for i in range(N_lam): 
    
            ## [Absorption and scattering for current wavelength]
            a = np.squeeze(optaa_a[:,lam_is[i]])
            b = np.squeeze(optaa_b[:,lam_is[i]])
         
            ## [Smoothing the absorption and scattering for current wavelength]
            z_s, a_s = ODF.Smooth_Profile_55(optaa_z, a)
            z_s, b_s = ODF.Smooth_Profile_55(optaa_z, b)
            
            ## [Chla being interpolated to absorption grid.]
            chla_s_interp = np.interp(z_s, z_chla_s, chla_s)
    
            ## [Loop over the phytoplankton species.]
            for k in range(N_phy): 
                ## [Filling the A array, recall we are solving for a single z location. 
                ##  Each row corresponds to a wavelength and each column a different species.]
                ## [Absorption array.]
                A[2*i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[0] * chla_s_interp[z_i])
                A[2*i+1,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[1] * chla_s_interp[z_i])
    
            ## [The rhs y, is simply the a_ooi minus the a_wat from Dut.]
            ## [The water absorption subtraction is commented out at the momment because it leads to negative.]

            ## [Get the cdom absorption]
            lam_cdom_i = ODF.Get_Wavelength_Index(wavelength_a[0,:], cdom_reflam)
            CDOM = OI.CDOM_refa(z_s, 
                                a_s[z_i], 
                                wavelength_a[0,lam_cdom_i],
                                wavelength_a[0,lam_is[i]])
            a_cdom = CDOM.a
    
            ## WARNING 
            ## WARNING The abs should be removed and the negative absorption problem should be fixed, this is temp. 
            ## WARNING 
            ## [The cdom is subtracted.]
            y[2*i] = abs(a_s[z_i]) - a_cdom # - abscat(wavelengths[i], 'water', C2chla='default')[0]
            y[2*i+1] = abs(b_s[z_i]) 

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
         
            ## [Smoothing the absorption and scattering for current wavelength]
            z_s, a_s = ODF.Smooth_Profile_55(optaa_z, a)
            z_s, b_s = ODF.Smooth_Profile_55(optaa_z, b)
    
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
            lam_cdom_i = ODF.Get_Wavelength_Index(wavelength_a[0,:], cdom_reflam)
            CDOM = OI.CDOM_refa(z_s, 
                                a_s[z_i], 
                                wavelength_a[0,lam_cdom_i],
                                wavelength_a[0,lam_is[i]])
            a_cdom = CDOM.a
    
            ## WARNING 
            ## WARNING The abs should be removed and the negative absorption problem should be fixed, this is temp. 
            ## WARNING 
            ## [The cdom is subtracted.]
            if  ab == 'a':
                y[i] = abs(a_s[z_i]) - a_cdom # - abscat(wavelengths[i], 'water', C2chla='default')[0]
            if  ab == 'b':
                y[i] = abs(b_s[z_i]) 

    ## [Solvving the problem as a convex optimization problem using the cvxpy.]
    ## [This is a least squares implementation fo the problem with constraints.]
    x = cp.Variable(N_phy)
    objective = cp.Minimize(cp.sum_squares(A@x - y))
    constraints = [0 <= x, cp.sum(x) == 1.0]
    prob = cp.Problem(objective, constraints)
    result = prob.solve()
    x = x.value

        
    ## [Solving the least square system for the phy ratios.]
    #x = np.linalg.lstsq(A, b)
#    x = scipy.optimize.nnls(A,y)[0]
#    x = scipy.optimize.lsq_linear(A,y, bounds=(0,1))

    ## [Divide x by its two norm.]
    ## [x should be positive.]
#    x = x/ sum(x)

    return A, y, x 


def Run_Lstsq_Time(N_profs, PI, wavelengths, z_i, phy_species, depth_profs, dt_profs, chla_profs, cdom_profs, 
                   optaa_dat, ab='ab'): 
    """
    Run the least square problem over many profiles in time. The output from this function 
    can then be used to construct a hoffmuller diagram of the species ratios.
    """

    ## [The number of phytoplankton species.]
    N_phy = len(phy_species)

    ## [The exmpty x_profs array, which contains the lst sq species estimate for each profile
    ##  at the same z_i.]
    x_profs = np.zeros((N_phy, N_profs))

    ## [The residuals of the fit for each.]
    residuals = np.zeros(N_profs)
    
    ## [Loop over the profiles.]
    for prof_index in range(N_profs):
        print(prof_index)
        ## [Running the least square estimation of the ratio of phytoplankton.]
        A, y, x = Est_Spec_Lstsq(PI, 
                                wavelengths,
                                z_i,
                                prof_index,
                                phy_species,
                                depth_profs, 
                                dt_profs, 
                                chla_profs, 
                                cdom_profs, 
                                optaa_dat, 
                                ab=ab,
                                plot=False)

        ## [If the point is bad such that y is all zeros.]
        if sum(np.abs(y)) == 0.0: 
            x = np.NaN * x

        ## [Storing the result into x_profs array.]
        x_profs[:,prof_index] = x

        ## [The two norm of the difference between y_true and the apporx solution.]
        residuals[prof_index] = np.linalg.norm(A@x - y, ord=2) / np.linalg.norm(y, ord=2)

    return x_profs, residuals


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

    if ab == 'a': 
        ab_str = 'Absorption'
    if ab == 'b': 
        ab_str = 'Scattering'

    if ab == 'ab':
        ## [Two subplots one for absorption and one for scattering.]
        fig, axs = plt.subplots(ncols=1, nrows=2)
        axs[1].set_xlabel("Wavelength [nm]")
        axs[0].set_title('Coupled Least Square Solution')

        ## [Loop over the scat/abs plot.]
        for k, ax in enumerate(axs): 

            ## [k==0 is the Absorption portion of the solution.]
            if k==0: 
                ab_str = 'Absorption'
            ## [k==0 is the Absorption portion of the solution.]
            if k==1: 
                ab_str = 'Scattering'

            ax.plot(wavelengths, y[k::2], 'k', label = f"OOI Total Effective {ab_str}")
            ax.plot(wavelengths, y_est[k::2], ':k', label = f"Least Square Approx {ab_str}")
         
            ## [Plot the species coefficients]
            for i in range(N_phy): 
                ax.plot(wavelengths, A[k::2, i], label = phy_species[i]) 
     
            ax.set_ylabel(f"Effective {ab_str}  [m^-1]")
     
            ax.grid()
            ax.legend()
        
        fig.show()

    ## [The uncoupled abs/scat specific solution.]
    else: 
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


def Plot_Spec_Lstsq_Hoffmuller(N_profs, dt_profs, phy_species, x_profs, plot_residuals=False, residuals=None): 
    """
    This plots a Hoffmuller diagram of the ratios of different species opf phytoplankton in time

    It will be done by plotting a stacked bar plot.]
    """

    ## [The number of species.]
    N_phy = len(phy_species)

    ## [Init the bottom coordinate for the stacked bar plot.]
    bot = np.zeros(N_profs)

    ## [Init the dt coordinate array.]
#    dt = np.zeros(N_profs, dtype=object)
    dt = np.zeros(N_profs)

    ## [The colors for each species.]
    colors = ['b','r','g','y','c','m','k','w'][:N_phy]

    ## [The time coordinate.]
    for k in range(N_profs):
#        dt[k] = dt_profs[k][0]
         dt[k] = k

    ## [Plot the residual on the same figure.]
    if plot_residuals: 
        fig, axs = plt.subplots(nrows=1, ncols=2)
    else:
        fig, ax = plt.subplots()

    ## [Top ax is the hoffmuller.]
    ax = axs[0]
    ## [Loop over the species.]
    for k in range(N_phy): 
        
        ## [plot the phy species layer of the bar plot.]
        ax.bar(dt, x_profs[k,:], bottom=bot, color = colors[k], label=phy_species[k], width=0.9)

        ## [Update the bottom of the bar plot.]
        bot = bot + x_profs[k,:]

    ax.set_title('Constrained Least Squares \n Phytoplankton Community Estimation')
    ax.set_ylabel('Fractional Concentraion')
    ax.set_xlabel('Temporal Index')
    ax.legend()
    ax.grid()
    fig.show()

    if plot_residuals: 
        ## [The bottom ax.]
        ax = axs[1]

        ## [Plot the residuals.]
        ax.bar(dt, residuals, color = 'k', width=0.9)

        ax.set_title('Two Norm of Residual Vector')
        ax.set_ylabel(r'$\frac{|\mathbf{A} \mathbf{x} - \mathbf{y} |_2}{|\mathbf{y}|_2}$')
        ax.set_xlabel('Temporal Index')
        ax.grid()
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
    ooi_savefile_head = "ooi_data/ooi_dat"

    ## [Data load flag, fale means load data from pickle.]
    download = False

    ## [Functional parameters.]
    N=100
    wavelengths = np.arange(425, 725, 25)
    #    phy_species = ['HLPro', 'LLPro', 'Cocco', 'Diat', 'Syn', 'Tricho', 'Lgeuk'] 
    phy_species = ['HLPro', 'Cocco', 'Diat', 'Syn'] 
    PI = Param_Init()
    cdom_reflam = 412.0
    prof_index = 1
    z_i = -2
    ## [The absorption or scattering flag.]
    ab = 'ab'
    ## [The number of profiles for the Hoffmuller.]
    N_profs = 50

    ## [Download or load the data sets.]
    ooi_data = ODF.Download_OOI_Data(ooi_savefile_head, 
                                     download, 
                                     site_name, 
                                     assembly, 
                                     method, 
                                     start, 
                                     stop)

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
                             z_i,
                             phy_species,
                             flort_prof, 
                             optaa_prof,
                             ab=ab,
                             plot=False)
    Plot_Abs_Est_LstSq_Sol(wavelengths, phy_species, A, y, x, ab=ab)

    ## [Running the least squares problem over many profiles.]
#    x_profs, residuals = Run_Lstsq_Time(N_profs, 
#                                        PI, 
#                                        wavelengths, 
#                                        z_i, 
#                                        phy_species, 
#                                        depth_profs, 
#                                        dt_profs, 
#                                        chla_profs, 
#                                        cdom_profs, 
#                                        optaa_dat, 
#                                        ab=ab)
#    Plot_Spec_Lstsq_Hoffmuller(N_profs, dt_profs, phy_species, x_profs, plot_residuals=True, residuals=residuals)
