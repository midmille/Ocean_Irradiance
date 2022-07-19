"""
Created: 29 November 2021
Author: Miles D Miller 

The implementation of the Baird et al. method for calculating the 
remote sensed reflectance given absorbtion and scattering profiles. 

The reference article is: 

Mark E. Baird 
27 November 2015
Remote-sensing reflectance and true colour produced by a coupled 
hydrodynamic, optical, sediment, biogeochemical model of the 
Great Barrier Reef, Australia: Comparison with satellite data

"""

## External
import numpy as np

## User
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat

def ocean_irradiance_baird(hbot, ab_wat, theta_air, phy = None, N = 30):

    """
    
    """

    Nm1 = N-1
    Nm2 = N-2

    ##unpacking the ab_wat_tuple 
    a_wat,b_wat = ab_wat 
    a = a_wat
    b = b_wat 
    
    ## If phytoplankton, otherwise just water in column.
    if phy: 
        
        ## unpacking the phytoplankton object
        z_phy = phy.z
        Nphy = phy.Nphy
        
        ## array for different phy
        phy_prof = phy.phy
        
        ## coefficients
        a_phy = phy.a
        b_phy = phy.b
        
        ## Just one phytoplankton species
        if Nphy == 1 : 
            a = a + phy_prof * a_phy
            b = b + phy_prof * b_phy
            
        ## More than one species
        elif Nphy > 1 : 
            for k in range(Nphy):
                a = a + phy_prof[:,k] * a_phy[k]  
                b = b + phy_prof[:,k] * b_phy[k]

    ## Log grid calculation. 
    ## The cell centers.
    z_e = OI.Log_Trans(hbot, N) 
    ## making the cell centers z. 
    z_c = np.zeros(Nm1)
    for k in range(Nm2, -1, -1): 
        z_c[k] = z_e[k+1] + (z_e[k+1] - z_e[k])/2 

    ## Interpolating a,b vectors from z_phy to z.
    if phy: 
        a = np.interp(z_c,z_phy,a)
        b = np.interp(z_c,z_phy,b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        
    ## Calculation of backscattering coefficient.
    b_b = .551*b
    b_f = b - b_b 
    
    #print('a', a)
    #print('b', b)


    ## Start of the Baird model

    ## The azimuth in water angle 
    theta_sw = np.arcsin((np.sin(theta_air))/1.33) 
    #print('theta_sw', theta_sw)

    ## vertical attenuation coefficient
    gi = 0.402
    gii = 0.180
    ## This K is on the cell centers.
    K_c = (a/(np.cos(theta_sw)))*(np.sqrt(1 + (gi + gii*(np.cos(theta_sw)))*(b/a)))
    ## This K is on the cell edges.
    K_e = np.interp(z_e, z_c, K_c)
    #print('K', K_e)

    ## The weighting function with the integral calculated using trap rule.
    ## Calculated on centers.
    w = np.zeros(Nm1)
    
    for k in range(Nm2, -1, -1):    
        w[k] = 0.5*(np.exp(-2*K_e[k+1]) + np.exp(-2*K_e[k])) 
    #print('w', w)

    ## Calculating u from equation 12
    u = 0
    for k in range(Nm2, -1, -1): 
     #   print(u)
        u = u + (((w[k]*b_b[k])/(a[k] + b_b[k])) * (z_e[k+1] - z_e[k]))

    ## rrs calculation.
    g0 = 0.0895
    g1 = 0.1247

    rrs = g0*u + g1*u**2
    #print(rrs)

    ## Rrs calculation. 
    Rrs = (0.52*rrs)/(1-1.7*rrs)

    return Rrs


if __name__ == '__main__': 

    wavelengths = [443, 551]

    ## N layers 
    N = 200

    ## azimuth angle in air, rn is just the overhead. 
    theta_air = .558
    ## The z grid for the phy prof
    z_phy = np.linspace(-600, 0, N)
    ## Artificial phytoplankton profile units chla. 
    phy_prof = OI.artificial_phy_prof(z_phy, 0, 10, 1)
    
    ## dictionary with lam as keys
    Rrs_dict = {}

    for lam in wavelengths:
        ab_wat = abscat(lam, 'water')
       
        a_phy, b_phy = abscat(lam, 'Diat')

        phy = OI.Phy(z_phy, phy_prof, a_phy, b_phy)

        hbot = z_phy[0] 

        Rrs_dict[lam] = ocean_irradiance_baird(hbot, ab_wat, theta_air, phy = phy, N = N)

    ## getting chla
    chla = OIR.OCx_alg(Rrs_dict[443], Rrs_dict[551]) 



 
