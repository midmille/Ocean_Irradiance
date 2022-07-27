"""
This file is to test the different ocean_irradiance methods against the analytical solution 
for different resolutions.
"""


## [User Modules.]
from ocean_irradiance_module import Ocean_Irradiance as OI
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
from ocean_irradiance_module.PARAMS import Param_Init
## [External Modules.]
import numpy as np
import matplotlib.pyplot as plt

def resolution_sensitivity(PI, methods, Ns, wavelengths):
    """

    
    """

    hbot = -500
    Nphy = 20
    z_phy = np.linspace(hbot, 0, Nphy)
    phy_type = 'Diat'
    phy_conc = 1
    phy_prof = np.full(Nphy,phy_conc) 

    ## [The linear grid]
    err_lin = np.zeros((len(wavelengths), len(methods), len(Ns)))
    ## [The log grid.]
    err_log = np.zeros((len(wavelengths), len(methods), len(Ns)))

    for i, lam in enumerate(wavelengths): 
        ## [The water absorption and scattering for lam.]
        ab_wat = abscat(lam, 'water')
        ab_phy = abscat(lam, phy_type)

        for k, method in enumerate(methods): 
            for j, N in enumerate(Ns): 
                phy = OI.Phy(z_phy, phy_prof, esd(phy_type), ab_phy[0], ab_phy[1])
                
                ## [The linear grid case.]
                PI.grid = 'linear'
                ## [Compute the analytical solution as the truth.]
                Eda, Esa, Eua, za = OI.ocean_irradiance(PI, hbot, ab_wat, method='analytical', phy=phy, N=N)
                ## [Compute the given method solution.]
                Ed, Es, Eu, z = OI.ocean_irradiance(PI, hbot, ab_wat, method=method, phy=phy, N=N)
                ## [Calculate the rmse.]
                err_lin[i, k, j] = np.linalg.norm((Esa - Es), ord=2) + np.linalg.norm((Eua-Eu), ord=2)

                ## [The log grid case.]
                PI.grid = 'log'
                ## [Compute the analytical solution as the truth.]
                Eda, Esa, Eua, za = OI.ocean_irradiance(PI, hbot, ab_wat, method='analytical', phy=phy, N=N)
                ## [Compute the given method solution.]
                Ed, Es, Eu, z = OI.ocean_irradiance(PI, hbot, ab_wat, method=method, phy=phy, N=N)
                ## [Calculate the rmse.]
                err_log[i, k, j] = np.linalg.norm((Esa - Es), ord=2) + np.linalg.norm((Eua-Eu), ord=2)

                

    return err_lin, err_log


def plot_resolution_sensitivity(PI, methods, Ns, wavelengths, err_lin, err_log): 
    """
    """
    
    fig, axs = plt.subplots(nrows = 1, ncols=len(wavelengths))

    for k, lam in enumerate(wavelengths): 
        ax = axs[k]

        for j, method in enumerate(methods): 
            
            l = ax.plot(Ns, err_lin[k,j,:], label=f'{method} linear grid')
            ax.plot(Ns, err_log[k,j,:], '--', label=f'{method} log grid', color = l.get_color())

        if k ==0: 
            ax.set_ylabel(r'Error [$|Es_{\mathrm{analytical}} - Es_{\mathrm{numerical}}|_2 + |Eu_{\mathrm{analytical}} - Eu_{\mathrm{numerical}}|_2$]')
            ax.legend()

        ax.set_title(f'Wavelength {lam} [nm]')
        ax.set_xlabel('Number of Grid Points')
        ax.grid()

    fig.show()


    return 


if __name__ == '__main__': 
    
    PI = Param_Init()
    PI.pt1_perc_zbot = False
    PI.pt1_perc_phy = False

    ## [The number of vertical layers used for each method.]
    Ns = np.arange(5, 20, 5)
    ## [The different methods]
    methods = ['shootdown', 'shootup', 'scipy', 'dutkiewicz']
    ## [The wavelengths for the study.]
    wavelengths = [443, 551]

    err_lin, err_log = resolution_sensitivity(PI, methods, Ns, wavelengths)

    plot_resolution_sensitivity(PI, methods, Ns, wavelengths, err_lin, err_log)