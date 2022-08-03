"""
This file is for the testing of the irradiance model's sensitivity to surface boundary conditions. Namely
the ratio of Ed0 to Es0.

"""

import numpy as np
import matplotlib.pyplot as plt

from ocean_irradiance_module import Ocean_Irradiance as OI 
from ocean_irradiance_module import Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module.PARAMS import Param_Init 
import ocean_irradiance_module.Wavelength_To_RGB as W2RGB

def Run_Study(Niv, Ed0s, Es0s, PI, hbot, N, wavelengths, z_phy, phy_prof, phy_type):
    """
    """

    Rrss_dict ={}
    for k, lam in enumerate(wavelengths):
        Rrss = np.zeros(Niv)
        ab_wat = abscat(lam, 'water')
        phy = OI.Phy(z_phy, phy_prof, esd(phy_type), abscat(lam, phy_type)[0], abscat(lam, phy_type)[1])
        for j in range(Niv): 
            print(j, '/', Niv)
            
            ## [Make the initital conditions the correct values.]
            PI.Ed0 = Ed0s[j]
            PI.Es0 = Es0s[j]

            Ed, Es, Eu, z = OI.ocean_irradiance(PI, 
                                             hbot, 
                                             ab_wat, 
                                             phy=phy, 
                                             N=N)
            ## [Calculate Rrs]
            Rrss[j] = OIR.R_RS(PI.Ed0, PI.Es0, Eu[-1])

        Rrss_dict[lam] = Rrss
        

    return Rrss_dict

def Plot_Result(wavelengths, Ed0s, Es0s, Rrss_dict): 
    """

    """

    fig, ax = plt.subplots()

    for k, lam in enumerate(wavelengths): 
        rgb = W2RGB.wavelength_to_rgb(lam)
        Rrss = Rrss_dict[lam]
        ax.plot(Ed0s, np.absolute(Rrss-Rrss[-1])/Rrss[-1], color=rgb, label =lam)

    ax.grid()
    ax.legend()

    fig.show()

        


    return 

if __name__ == '__main__':

    PI = Param_Init()

    PI.pts_perc_phy = True
    PI.pts_perc_zbot = True
    
    PI.grid = 'log'
    
    ##[The resolution.]
    N = 100

    wavelengths = PI.wavelengths

    hbot = -200

    phy_type = 'Diat'
    z_phy = np.linspace(hbot,0,N)
    phy_conc = 10
    phy_prof = np.full(N, phy_conc)
    
    ## [The number of InitiL VALUES TO TRY.]
    Niv = 10

    Ed0s = np.linspace(0, 1, Niv)
    Es0s = 1-Ed0s

    Rrss_dict  = Run_Study(Niv, Ed0s, Es0s, PI, hbot, N, wavelengths, z_phy, phy_prof, phy_type)

    Plot_Result(wavelengths, Ed0s, Es0s, Rrss_dict)
