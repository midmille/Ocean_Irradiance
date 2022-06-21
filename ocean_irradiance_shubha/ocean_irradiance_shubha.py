"""
This is to implement the two stream model provided by shubha 
The translation of her model to our variables was made on February 16 2022
"""
## external Mods
import numpy as np
from netCDF4 import Dataset 
## user mods
from ocean_irradiance_module.PARAMS import Param_Init 
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as ESD
from ocean_irradiance_module import Ocean_Irradiance as OI
#from ocean_irradiance_module import Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module import Wavelength_To_RGB
#import os
#import sys
#import pickle
#import matplotlib as mpl 
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature



def numerical_Ed(z, a, b_b, v_d, Ed0):
    
    ## The coefficient of attentuation
    c_d = (a+b_b) / v_d
    N = len(z)
    Nm1 = N - 1
    
    Ed = np.zeros(N)
    Ed[Nm1] = Ed0

    for k in range(Nm1-1, -1, -1) :
        
        dz = z[k] - z[k+1]
        dzo2 = dz / 2
        
        dEddz1 = c_d[k]*Ed[k+1]
        dEddz2 = c_d[k]*(Ed[k+1]+(dEddz1*dzo2))
        dEddz3 = c_d[k]*(Ed[k+1]+(dEddz2*dzo2))
        dEddz4 = c_d[k]*(Ed[k+1]+(dEddz3*dz))
        Ed[k] = Ed[k+1] + (( (dEddz1/6)+(dEddz2/3)+(dEddz3/3)+(dEddz4/6) )*dz)
        
    return Ed


def numerical_Eu(z, Ed, a, b_b, v_u, v_d): 

    ## calculating Eu from shubhas two stream model.  
    N = len(z)
    Nm1 = N - 1
    
    Eu = np.zeros(N)
    ## Assuming Eu at z=0 is also 0
    Eu[0] = 0
    
    ## Looping water column 
    for k in range(0, Nm1): 
        dz = z[k] - z[k+1]
        dzo2 = dz / 2
        
        dEudz1 = ((a[k]+b_b[k])/v_u)*Eu[k] - (b_b[k]/v_d)*Ed[k]
        dEudz2 = ((a[k]+b_b[k])/v_u)*(Eu[k]+(dEudz1*dzo2)) - (b_b[k]/v_d)*Ed[k]
        dEudz3 = ((a[k]+b_b[k])/v_u)*(Eu[k]+(dEudz2*dzo2)) - (b_b[k]/v_d)*Ed[k]
        dEudz4 = ((a[k]+b_b[k])/v_u)*(Eu[k]+(dEudz2*dz)) - (b_b[k]/v_d)*Ed[k]
        Eu[k+1] = Eu[k] + (( (dEudz1/6)+(dEudz2/3)+(dEudz3/3)+(dEudz4/6) )*dz)

    return Eu


def ocean_irradiance_two_stream_ab(hbot, ab_wat, N, phy=None, CDOM_sal=None, CDOM_dens=None, CDOM_refa=None, pt1_perc_zbot=True, pt1_perc_phy=True):

    """
    This function calculates the irradiance grid, calculates the absorption, scattering, and backscatter 
    profiles and then interpolates those onto the irradiance grid

    Parameters
    ----------

    Returns
    -------
    """

    ##unpacking the ab_wat_tuple 
    a_wat,b_wat = ab_wat 
    a = a_wat
    b = b_wat 
    
    ## backscattering due to phy
    b_b_phy = 0
    ## If phytoplankton, otherwise just water in column.
    if phy: 
        print('HERE')
        ## unpacking the phytoplankton object
        z_phy = phy.z
        Nphy = phy.Nphy
        
        ## array for different phy
        phy_prof = phy.phy
        
        print('phy_prof:', phy_prof)
        ## coefficients
        a_phy = phy.a
        b_phy = phy.b
        
        ## equivalent spherical diameter
        esd = phy.esd
        ## Just one phytoplankton species
        if Nphy == 1 : 
            ## The back scatter ratio
            bb_r = OI.Backscatter_Ratio(esd)    
            print('bb_r', bb_r)
            a = a + phy_prof * a_phy
            b = b + phy_prof * b_phy
            b_b_phy = b_b_phy + phy_prof * b_phy * bb_r
            
        ## More than one species
        elif Nphy > 1 : 
            for k in range(Nphy):
                ## The back scatter ratio
                bb_r = Backscatter_Ratio(esd[k])    
                a = a + phy_prof[:,k] * a_phy[k]  
                b = b + phy_prof[:,k] * b_phy[k]
                b_b_phy = b_b_phy + phy_prof[:,k] * b_phy[k] * bb_r

    ## Inclusion of CDOM
    if CDOM_sal: 
        ## unpacking the object
        ## For now it is assumed that phy and cdom share same grid.
        if phy:
            a_cdom = CDOM_sal.a
            a = a + a_cdom 
        ## This is for only CDOM no phy
        else:
            a_cdom = CDOM_sal.a
            z_cdom = CDOM_sal.z 
            a = a + a_cdom 

    ## [Inclusion of CDOM via its observed concentration.]
    if CDOM_dens: 
        ## [Unpack the object.]
        z_cdom = CDOM_dens.z 
        cdom = CDOM_dens.cdom
        cdom_a =  CDOM_dens.a 
        CDOM2C = CDOM_dens.CDOM2C

        ## [Interpolate CDOM to the phytoplankton grid.]
        if phy: 
            ## [Now cdom is on the phy grid.]
            a_cdom = cdom_a * CDOM2C * np.interp(z_phy, z_cdom, cdom)
            a += a_cdom
        ## [If there is not phy.]
        else: 
            ## [On the cdom grid.]
            a_cdom = cdom_a * CDOM2C * cdom 
            a += a_cdom

    ## [Inclusion of CDOM via its reference absorption.]
    if CDOM_refa: 
        ## [Unpack.]
        z_cdom = CDOM_refa.z
        a_cdom = CDOM_refa.a
        ## [Interpolate CDOM to the phytoplankton grid.]
        if phy: 
            ## [Now cdom is on the phy grid.]
            a_cdom = np.interp(z_phy, z_cdom, a_cdom)
            a += a_cdom
        ## [If there is not phy.]
        else: 
            ## [On the cdom grid.]
            a += a_cdom

    ## [One estimation of CDOM at a time.]
    if CDOM_sal and CDOM_dens:
        raise Exception("Can't have both CDOM versions at same time")
    if CDOM_sal and CDOM_refa:
        raise Exception("Can't have both CDOM versions at same time")
    if CDOM_dens and CDOM_refa:
        raise Exception("Can't have both CDOM versions at same time")

    
    ## Irradiance Grid Stuff
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        if pt1_perc_phy == True:
            c_d = (a+b)/v_d
            zbot_pt1perc = OI.zbot_func(Ed0, a, b, v_d, phy=True, z=z_phy) 
        else:    
            c_wat = (a_wat + b_wat)/v_d
            zbot_pt1perc = OI.zbot_func(Ed0, a_wat, b_wat, v_d, phy=False)
        if zbot_pt1perc == None:
            print('bad pt1 perc light level')
            zbot_pt1perc = -100
        #print(zbot_pt1perc)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    z = OI.Log_Trans(zbot, N) 
    ## linear z 
    #z = np.linspace(zbot, 0, N)
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    if phy: 
        print('a_interp_pre:', a)
        print('z-phy_interp:', z_phy)
        print('z_interp:', z)
        a = np.interp(z,z_phy,a)
        print('a_interp_post:', a)
        b = np.interp(z,z_phy,b)
        b_b_phy = np.interp(z, z_phy, b_b_phy)
    elif CDOM_sal or CDOM_dens or CDOM_refa: 
        a = np.interp(z,z_cdom,a)
        ## no scattering for cdom, thus just water scattering.
        b = np.full(N, b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        
    b_b_wat = .551*b_wat

    b_b = b_b_wat + b_b_phy

    return z, a, b, b_b
 
def ocean_irradiance_shubha(hbot, Ed0, ab_wat, coefficients, phy=None, CDOM=None, N=30, 
                            pt1_perc_zbot=True, pt1_perc_phy=True, a=None, b_b=None):

    """
    The implementation of the two stream model of Shubha 1997, 1998. 
    """

    ## PARAMS FROM DUTKIEWICZ 2015 
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    ##N centers
    Nm1 = N - 1  

    Ed1 = np.zeros(N)
    Eu1 = np.zeros(N) 
    
    ## BC
    Ed1[Nm1] = Ed0 ##surface 
    Eu1[0] = 0 ##bottom
   
    Ed=np.copy(Ed1)
    Eu=np.copy(Eu1)

    ## [Calculate the irradiance grid, a, b, b_b.]
    z, a, b, b_b = ocean_irradiance_two_stream_ab(hbot, ab_wat, N, phy=phy, CDOM=CDOM, 
                                                  pt1_perc_zbot=pt1_perc_zbot,
                                                  pt1_perc_phy=pt1_perc_phy)


    Ed = numerical_Ed(z, a, b_b, v_d, Ed0)

    Eu = numerical_Eu(z, Ed, a, b_b, v_u, v_d)

    return Ed, Eu, z, a, b_b


def ocean_irradiance_shubha_ooi(hbot, Ed0, coefficients, z_a, a, z_b_b, b_b, N=30, pt1_perc_zbot=True):

    """
    The implementation of the two stream model of Shubha 1997, 1998. 
    
    This implementation is different from the function above mainly because it does not calculate the 
    absorption and scattering on its own. It has the absorption and scattering arrays as input. 
    
    This was created to be used with the ooi data set which gives the absorption and scattering arrays. 
    """

    ## PARAMS FROM DUTKIEWICZ 2015 
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    ##N centers
    Nm1 = N - 1  

    Ed1 = np.zeros(N)
    Eu1 = np.zeros(N) 
    
    ## BC
    Ed1[Nm1] = Ed0 ##surface 
    Eu1[0] = 0 ##bottom
    
    ## Irradiance Grid Stuff
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        c_d = (a+b)/v_d
        zbot_pt1perc = OI.zbot_func(Ed0, a, b, v_d, phy=True, z=z_phy) 
        if zbot_pt1perc == None:
            print('bad pt1 perc light level')
            zbot_pt1perc = -100
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 

    ## log transformed z grid.
    z = OI.Log_Trans(zbot, N) 
    ## linear z 
    #z = np.linspace(zbot, 0, N)
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    a = np.interp(z,z_a,a)
    b_b = np.interp(z,z_b_b,b_b)
        
    Ed=np.copy(Ed1)
    Eu=np.copy(Eu1)

    Ed = numerical_Ed(z, a, b_b, v_d, Ed0)

    Eu = numerical_Eu(z, Ed, a, b_b, v_u, v_d)

    return Ed, Eu, z
           
        
def Demo(): 
    ## Demonstration of finding the irradiance profiles for an artificial 
    ## phytoplankton concentration profile
    
    ## Absorbtion and scattering coefficients are taken from Dutkiwicz 2015 
    ## for the Diatom species. 
    import matplotlib.pyplot as plt
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    PI = Param_Init()
    
    N = 200
    Nm1 = N-1 
    lams = [443, 551]
    
    z = np.linspace(-200,0,N)

    phy_prof = OI.artificial_phy_prof(z, -10, 20, 1, prof_type = 'gauss')
    #phy_prof = np.full(200, 1) 
    #phy_prof = np.append(phy_prof1, np.full(100, 1))
    # ROMS_point = np.genfromtxt('ChrisData_good_point.csv', delimiter=',')
    # phy_prof = ROMS_point[1:,2]
    # print(phy_prof)
    # z = ROMS_point[1:,0]
    # phy_prof = np.full(len(z), 1)
    fig, axes = plt.subplots(1, 3, sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    
    phy_type = 'Syn'
    for k, lam in enumerate(lams):
        ab_wat = abscat(lam, 'water')
    
        a_phy, b_phy = abscat(lam, phy_type)
        
        
        ## Define the Phytoplankton class.
        phy = OI.Phy(z, phy_prof, ESD(phy_type),a_phy, b_phy)
    
        ## Salinity for CDOM
        salt = np.full(N, 34.5)
        ## define the CDOM class
        cdom = OI.CDOM(z,salt,lam) 
        
    
        ## The fixed point position: 
        fp = -50
        fpi =0
    
        zbot = z[0]

        Ed, Eu, zarr = ocean_irradiance_shubha(zbot, 1, ab_wat, PI.coefficients, phy=phy, CDOM=None, N=N, pt1_perc_zbot=True, pt1_perc_phy=False)
       
        ## Dutkiewicz model.
        Ed_d, Es_d, Eu_d, z_d = OI.ocean_irradiance_shoot_up(zbot,1,0,PI.Euh,ab_wat, PI.coefficients, phy=phy, CDOM=None, N=N, pt1_perc_zbot = False, pt1_perc_phy=False)
        print(Es_d)


            
         
        ## Plotting the Results
        #-------------------------------------------------------------------------
        
        markers = ['-', '-'] 
        Ed_c = 'g'
        Es_c = 'b'
        Eu_c = 'r'
        if lam == 443:
            ax2.plot(Ed, zarr, label=f'Ed', color = Ed_c, ls=markers[k] )
            ax2.plot(Eu, zarr, label=f'Eu', color = Eu_c, ls= markers[k])
            ## plotting the dtukiewicz model also. 
#            ax2.plot(Ed_d + Es_d, z_d, label=f'Ed+Es Dut', color = Ed_c, ls=':' )
#            ax2.plot(Eu_d, z_d, label=f'Eu Dut', color = Eu_c, ls=':' )
            
            ## calculate chl-a
            #Rrs443 = OIR.R_RS(PI.Ed0, PI.Es0, Eu[-1])
        if lam == 551:
            ax3.plot(Ed, zarr, label=f'Ed', color = Ed_c, ls=markers[k] )
            ax3.plot(Eu, zarr, label=f'Eu', color = Eu_c, ls= markers[k])
            ## plotting the dtukiewicz model also. 
#            ax3.plot(Ed_d +Es_d, z_d, label=f'Ed+Es Dut', color = Ed_c, ls=':' )
#            ax3.plot(Eu_d, z_d, label=f'Eu Dut', color = Eu_c, ls=':' )

            #Rrs551 = OIR.R_RS(PI.Ed0, PI.Es0, Eu[-1])

    ax1.plot(phy_prof, z)
    ax1.set_xlabel('Concentration [mg Chl-a m^-3]')
    ax1.set_ylabel('Z [m]')
    ax1.set_title('Phytoplankton Concentration')
    ax1.grid()
 
    ax2.set_xlim(-0.1,1)
    ax2.set_xlabel('Irradiance')
    ax2.set_title(r'Irradiance Profiles $\lambda=443$')
    ax2.legend()
    ax2.grid()

    ax3.set_xlim(-0.1,1)
    ax3.set_xlabel('Irradiance')
    ax3.set_title(r'Irradiance Profiles $\lambda=551$')
    ax3.legend()
    ax3.grid()
    
    fig.show()
    
    ## caclulate chla
    #chla = OIR.OCx_alg(Rrs443, Rrs551)
    print('chla', chla)
    
    return zarr, Ed, Eu


if __name__ == '__main__': 


    Demo()










