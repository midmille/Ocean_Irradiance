"""
This file is for plotting the spectral absorption and scattering. It is meant to mimic figure one
in Dutkiewicz et al (2015)
"""

## [User Modules.]
from ocean_irradiance_module import Ocean_Irradiance as OI
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.Phytoplankton_Colormap import Get_Phy_Cmap_Dict

## [External Modules.]
import numpy as np
import matplotlib.pyplot as plt


def plot_constituent_abscat(phy_species): 
    """
    """
    ## [The number of wavelenghts.]
    Nlam = 100

    ## [The wavelengths for the plot.]
    wavelengths = np.linspace(400, 700, Nlam)

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize = (16, 6))

    ## [The water absorption and scattering.]
    a_wat = np.zeros(Nlam)
    b_wat = np.zeros(Nlam)

    for k, lam in enumerate(wavelengths): 
        a_wat[k] = abscat(lam, 'water')[0]
        b_wat[k] = abscat(lam, 'water')[1]

    ax_awat = axs[0]
    ax_bwat = ax_awat.twinx()

    ln0 = ax_awat.plot(wavelengths, a_wat, color = 'b', label = 'Absorption')
    ln1 = ax_bwat.plot(wavelengths, b_wat, color = 'g', label = 'Scattering')

    ax_awat.set_ylabel(r'Absorption [$\mathrm{m}^{-1}$]')
    ax_bwat.set_ylabel(r'Scattering [$\mathrm{m}^{-1}$]')
    ax_awat.set_xlabel('Wavelength [nm]')
    ax_awat.set_title('a) Water Absorption and Scattering')
    lns = ln0+ln1
    labs = [l.get_label() for l in lns]
    ax_awat.legend(lns, labs)
    ax_awat.grid()

    
    ## [The phytoplankton absorption and scattering.]
    ax_aphy = axs[1]
    ax_bphy = axs[2]

    ## [The color map dict for the different species of phytoplankton.] 
    cmap = Get_Phy_Cmap_Dict()

    ## [Loop the species.]
    for k, phy_type in enumerate(phy_species): 
        
        a_phy = np.zeros(Nlam)
        b_phy = np.zeros(Nlam)

        for j, lam in enumerate(wavelengths): 
            a_phy[j] = abscat(lam, phy_type)[0]    
            b_phy[j] = abscat(lam, phy_type)[1]    

        ax_aphy.plot(wavelengths, a_phy, color=cmap[phy_type], label=phy_type)
        ax_bphy.plot(wavelengths, b_phy, color=cmap[phy_type], label=phy_type)
    
    ## [labels for phytoplankton plots.]
    ax_aphy.set_ylabel(r'Absorption [$\mathrm{m}^2 \mathrm{mgChl}^{-1}$]')
    ax_aphy.set_xlabel('Wavelength [nm]')
    ax_aphy.set_title('b) Phytoplankton Specific Absorption')
    ax_aphy.legend()
    ax_aphy.grid()

    ax_bphy.set_ylabel(r'Scattering [$\mathrm{m}^2 \mathrm{mgChl}^{-1}$]')
    ax_bphy.set_xlabel('Wavelength [nm]')
    ax_bphy.set_title('c) Phytoplankton Specific Scattering')
    #ax_bphy.legend()
    ax_bphy.grid()
    
    plt.tight_layout()

    fig.show()
        
    return


if __name__ == '__main__': 

    PI = Param_Init()

    phy_species = PI.phy_species

    plot_constituent_abscat(phy_species)
