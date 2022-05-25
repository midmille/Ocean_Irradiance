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


def Est_Spec_Lstsq(PI, wavelengths, depthz, phy_species, flort_prof, optaa_prof, ab='a', plot=False,
                   weight=None): 
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

    ## [Get the depth index for optaa.]
    zi_optaa = (optaa_prof['depth'] - depthz).argmin()
    ## [Get the depth index for flort.]
    zi_flort = (flort_prof['depth'] - depthz).argmin()

    ## [The reference wavelength for the cdom absorption.]
    cdom_reflam = 412

    ## [This retreives the profiles]
    optaa_dt = optaa_prof['time'].data
    optaa_z = -optaa_prof['depth'].data
    ## !!!!!!!!!!!!!!!!TRANSPOSE IS TEMPORARY PUT THIS EARLIER ON IN OOI_DATA_FUNCTIONS
    ## !!!!!!!!!!!!!!!!TRANSPOSE IS TEMPORARY PUT THIS EARLIER ON IN OOI_DATA_FUNCTIONS
    ## !!!!!!!!!!!!!!!!TRANSPOSE IS TEMPORARY PUT THIS EARLIER ON IN OOI_DATA_FUNCTIONS
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
    flort_z = -flort_prof['depth'].data
    flort_chla = flort_prof['fluorometric_chlorophyll_a'].data
    z_chla_s, chla_s = ODF.Smooth_Profile_55(flort_z, flort_chla)

    ## [The wieghting for the leats squares problem.]
    ## [If the weights are provided are argument, unpack and use.]
    if weight: 
        weighta, weightb = weight 
    ## [The default weight is the ratio of spectral mean absorption to spectral mean scattering.]
    else:
        ## [Since the weiht is a ratio of a/b, then only multiply scattering by weight.]
        weighta = 1 / np.nanmean(optaa_a[zi_optaa,:])
        weightb = 1 / np.nanmean(optaa_b[zi_optaa,:])

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
            a_s = a
            b_s = b 
            z_s = optaa_z
         
            ## [Smoothing the absorption and scattering for current wavelength]
#            z_s, a_s = ODF.Smooth_Profile_55(optaa_z, a)
#            z_s, b_s = ODF.Smooth_Profile_55(optaa_z, b)
            
            ## [Chla being interpolated to absorption grid.]
            chla_s_interp = np.interp(z_s, z_chla_s, chla_s)

            ## [Get the cdom absorption]
            lam_cdom_i = ODF.Get_Wavelength_Index(wavelength_a[0,:], cdom_reflam)
            CDOM = OI.CDOM_refa(z_s, 
                                a_s[zi_optaa], 
                                wavelength_a[0,lam_cdom_i],
                                wavelength_a[0,lam_is[i]])
            a_cdom = CDOM.a
    
            ## [Loop over the phytoplankton species.]
            for k in range(N_phy): 
                ## [Filling the A array, recall we are solving for a single z location. 
                ##  Each row corresponds to a wavelength and each column a different species.]
                ## [Absorption array.]
                ## [Multiply the absorption eqaution by weight a that is (1/(avg(absorptio(lam)))).]
                ## [Multiply the scattering equation by weight b. Which is (1/avg(scattering(lam))).] 
                A[2*i,k] = ((abscat(wavelengths[i], phy_species[k], C2chla='default')[0] * chla_s_interp[zi_optaa])) * weighta 

                A[2*i+1,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[1] * chla_s_interp[zi_optaa]) * weightb
    
            ## [The rhs y, is simply the a_ooi minus the a_wat from Dut.]
            ## [The water absorption subtraction is commented out at the momment because it leads to negative.]
    
            ## WARNING 
            ## WARNING The abs should be removed and the negative absorption problem should be fixed, this is temp. 
            ## WARNING 
            ## [The cdom is subtracted.]
            y[2*i] = (abs(a_s[zi_optaa]) - a_cdom) * weighta  
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
                    A[i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[0] * chla_s_interp[zi_flort])
                if ab == 'b':
                    A[i,k] = (abscat(wavelengths[i], phy_species[k], C2chla='default')[1] * chla_s_interp[zi_flort])
    
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

    ## [Solvving the problem as a convex optimization problem using the cvxpy.]
    ## [This is a least squares implementation fo the problem with constraints.]
    x = cp.Variable(N_phy)
    objective = cp.Minimize(cp.sum_squares(A@x - y))
    constraints = [0 <= x, cp.sum(x) == 1.0]
    prob = cp.Problem(objective, constraints)
    result = prob.solve()
    x = x.value

    ## [If coupled then remove weight from returned array.]
    if ab == 'ab':
        A[0:2*N_lam:2, :] = A[0:2*N_lam:2, :] / weighta
        A[1:2*N_lam:2, :] = A[1:2*N_lam:2, :] / weightb
        y[0:2*N_lam:2] = y[0:2*N_lam:2] / weighta
        y[1:2*N_lam:2] = y[1:2*N_lam:2] / weightb
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

    ## [Text box params for the plot.]
    props = dict(facecolor='grey', alpha=0.6)
    txt = [f'{phy_species[k]}: {round(x[k],2)}' for k in range(len(phy_species))]
    txt = tuple(txt)
    txt = '\n'.join(txt)

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
                ax.plot(wavelengths, A[k::2, i], label = f'{phy_species[i]}: {round(x[i],2)}') 
     
            ax.set_ylabel(f"Effective {ab_str}  [m^-1]")
     
            ax.grid()
            ax.legend()

        ## [Setting som etext with the reslting coefficients.]
#        ylim = axs[0].get_ylim()
#        xlim = axs[0].get_xlim()
#        axs[0].text(xlim[0] + (xlim[1]-xlim[0])*0.1, ylim[0] + (ylim[1]-ylim[0])*0.1, txt, bbox=props)
        
        fig.show()

    ## [The uncoupled abs/scat specific solution.]
    else: 
        fig, ax = plt.subplots()
        ax.plot(wavelengths, y, 'k', label = f"OOI Total Effective {ab_str}")
        ax.plot(wavelengths, y_est, ':k', label = f"Least Square Approx {ab_str}")
     
        ## [Plot the species coefficients]
        for k in range(N_phy): 
            ax.plot(wavelengths, A[:, k], label = f'{phy_species[k]}: {round(x[k],2)}') 
    
        ax.set_ylabel(f"Effective {ab_str}  [m^-1]")
        ax.set_xlabel("Wavelength [nm]")
    
        ax.grid()
        ax.legend()

        ## [Setting som etext with the reslting coefficients.]
#        ylim = ax.get_ylim()
#        xlim = ax.get_xlim()
#        ax.text(xlim*0.1, ylim*0.1, txt, bbox=props)
    
        fig.show()
    
    return 


def Plot_Spec_Lstsq_Hoffmuller(N_profs, depthz, flort_profs, optaa_profs, phy_species, x_profs, plot_residuals=False, residuals=None): 
    """
    This plots a Hoffmuller diagram of the ratios of different species opf phytoplankton in time

    It will be done by plotting a stacked bar plot.]
    """

    ## [The number of species.]
    N_phy = len(phy_species)


    ## [Init the bottom coordinate for the stacked bar plot.]
    bot = np.zeros(N_profs)

    ## [The colors for each species.]
    colors = ['b','r','g','y','c','m','k','w'][:N_phy]

    ## [Init the dt coordinate array.]
    dt = []
    dt_ticks = []
    dt_labs = []
    ## [The label font size.]
    xlfs = 7
    ## [The time coordinate and their accompanying labels.]
    same_date = 0
    for k in range(N_profs):
        dt.append(k)
        ## [If the date is different then previously then set the date as the ticker.]
        if (str(optaa_profs[k]['time'].data.astype('datetime64[D]')[0])[-2:] != same_date): 
            dt_labs.append(str(optaa_profs[k]['time'].data.astype("datetime64[D]")[0]))
            dt_ticks.append(k)
            same_date = str(optaa_profs[k]['time'].data.astype('datetime64[D]')[0])[-2:]
        ## [Else use the time as the ticker.]
        else: 
            dt_labs.append(str(optaa_profs[k]['time'].data.astype("datetime64[m]")[0])[-5:])
            dt_ticks.append(k)

    ## [Plot the residual on the same figure.]
    if plot_residuals: 
        fig, axs = plt.subplots(nrows=1, ncols=2)
        ## [Least square axis.]
        ax_ls = axs[0]
        ## [Chla twin axis on least squares plot.]
        ax_cl = ax_ls.twinx()
        ## [Residual axis.]
        ax_rs = axs[1]
    else:
        fig, ax = plt.subplots()
        ## [Least square axis.]
        ax_ls = ax
        ## [Chla twin axis on least squares plot.]
        ax_cl = ax_ls.twinx()
 
    ## [Loop over the species.]
    for k in range(N_phy): 
        
        ## [plot the phy species layer of the bar plot.]
        ax_ls.bar(dt, x_profs[k,:], bottom=bot, color = colors[k], label=phy_species[k], width=0.9)

        ## [Update the bottom of the bar plot.]
        bot = bot + x_profs[k,:]

    ax_ls.set_title(f'Constrained Least Squares \n Phytoplankton Community Estimation at Depth: {depthz}[m]')
    ax_ls.set_ylabel('Fractional Concentraion')
    ax_ls.set_xlabel('Date Time')
    ax_ls.set_xticks(dt_ticks)
    ax_ls.set_xticklabels(dt_labs, rotation=75, fontsize=xlfs)
    ax_ls.legend(loc=2)
    ax_ls.grid(axis='y')

    ## [Getting the chla at the givenz level in each profile.]
    chla = np.zeros(N_profs)
    for k in range(N_profs):
        ## [Get the depth index for flort.]
        zi_flort = (flort_profs[k]['depth'] - depthz).argmin()
        chla[k] = flort_profs[k]['fluorometric_chlorophyll_a'].data[zi_flort]

    ## [Adding chla as twin plot to least squares.]
    ax_cl.plot(dt, chla, 'k', label='Chl-a', linewidth=2)
    ax_cl.set_yticks(np.linspace(ax_cl.get_yticks()[0], ax_cl.get_yticks()[-1], len(ax_ls.get_yticks())))
    ax_cl.set_ylabel(f"Fluorometirc Chl-a [{flort_profs[0]['fluorometric_chlorophyll_a'].attrs['units'][0]}]")
    ax_cl.legend(loc=1)

    if plot_residuals: 
        ## [Plot the residuals.]
        ax_rs.bar(dt, residuals, color = 'k', width=0.9)

        ax_rs.set_title('Two Norm of Residual Vector')
        ax_rs.set_ylabel(r'$\frac{|\mathbf{A} \mathbf{x} - \mathbf{y} |_2}{|\mathbf{y}|_2}$')
        ax_rs.set_xlabel('Date Time')
        ax_rs.set_xticks(dt_ticks)
        ax_rs.set_xticklabels(dt_labs, rotation=75, fontsize=xlfs)
        ax_rs.set_xticklabels(dt_labs)
        ax_rs.grid(axis='y')

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
    optaa_profs_save = "ooi_data/ooi_dat_rawoptaaprofs.p"

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
    depthz = 5
    ## [The absorption or scattering flag.]
    ab = 'ab'
    ## [The number of profiles for the Hoffmuller.]
    N_profs = 30

    ## [Download or load the data sets.]
    ooi_data = ODF.Download_OOI_Data(ooi_savefile_head, 
                                     download, 
                                     site_name, 
                                     assembly, 
                                     method, 
                                     start, 
                                     stop)

    flort_dat, spkir_dat, optaa_dat, flort_profs, spkir_profs, optaa_profs = ooi_data

    ## [Use the processed raw optaa profs.]
    ## !!!!! ALSO TEMPORARY PUT EARLIER ON...
    ## !!!!! ALSO TEMPORARY PUT EARLIER ON...
    ## !!!!! ALSO TEMPORARY PUT EARLIER ON...
    optaa_profs = pickle.load(open(optaa_profs_save,'rb'))
                                                            
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
                             ab=ab,
                             plot=False, 
                             weight = [0.4, 0.04])

    Plot_Abs_Est_LstSq_Sol(wavelengths, phy_species, A, y, x, ab=ab)

    ## [Running the least squares problem over many profiles.]
    x_profs, residuals = Run_Lstsq_Time(N_profs, 
                                        PI, 
                                        wavelengths, 
                                        depthz, 
                                        phy_species, 
                                        flort_profs, 
                                        optaa_profs,
                                        ab=ab)
    Plot_Spec_Lstsq_Hoffmuller(N_profs, depthz, flort_profs, optaa_profs, phy_species, x_profs, plot_residuals=True, residuals=residuals)
