"""
plot_est_ab_ooi.py
Author : Miles D. Miller, Univeristy of California Santa Cruz
Created : May 27, 2022 12:13pm
About : This file is for the plotting of the results of the functions and their 
        subsequent implementations in est_ab_ooi.py. Plotting the resulting least
        squares estimations and the Hoffmuller Diagram.   
"""

## [External Modules.]
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def Plot_Abs_Est_LstSq_Sol(wavelengths, optaa_prof, phy_species, A, y, x, phycolors, ab='a'): 
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

    ## [The profile start date for the title.]
    prof_sd = str(optaa_prof['time'].data.astype('datetime64[D]')[0])
    prof_t = str(optaa_prof['time'].data.astype('datetime64[m]')[0])[-5:]
 
    if ab == 'ab':
        ## [Two subplots one for absorption and one for scattering.]
        fig, axs = plt.subplots(ncols=1, nrows=2)
        axs[1].set_xlabel("Wavelength [nm]")
        axs[0].set_title(f'Constrained Least Squares Phytoplankton \n Community Estimation at Profile Date: {prof_sd}, Time: {prof_t}')
 
        ## [Loop over the scat/abs plot.]
        for k, ax in enumerate(axs): 

            ## [k==0 is the Absorption portion of the solution.]
            if k==0: 
                ab_str = 'Absorption'
            ## [k==0 is the Absorption portion of the solution.]
            if k==1: 
                ab_str = 'Scattering'

            ## [Plot the species coefficients]
            for i in range(N_phy): 
                ax.plot(wavelengths, A[k::2, i], lw='2', color = phycolors[phy_species[i]], label = f'{phy_species[i]}: {abs(round(x[i],2))}') 

            ax.plot(wavelengths, y[k::2], 'k', lw='2', label = f"OOI Total Effective {ab_str}")
            ax.plot(wavelengths, y_est[k::2], 'k', lw='2', linestyle='dotted', label = f"Least Square Approx {ab_str}")
     
            ax.set_ylabel(f"Effective {ab_str}  [m^-1]")
     
            ax.grid()
            if k==0:
                ax.legend()

        ## [Setting som etext with the reslting coefficients.]
#        ylim = axs[0].get_ylim()
#        xlim = axs[0].get_xlim()
#        axs[0].text(xlim[0] + (xlim[1]-xlim[0])*0.1, ylim[0] + (ylim[1]-ylim[0])*0.1, txt, bbox=props)
        
        fig.show()

    ## [The uncoupled abs/scat specific solution.]
    else: 
        fig, ax = plt.subplots()
     
        ## [Plot the species coefficients]
        for k in range(N_phy): 
            ax.plot(wavelengths, A[:, k], lw='2', color= phycolors[phy_species[k]], label = f'{phy_species[k]}: {abs(round(x[k],2))}') 
    
        ax.plot(wavelengths, y, 'k', lw='2', label = f"OOI Total Effective {ab_str}")
        ax.plot(wavelengths, y_est, 'k', lw='2', linestyle='dotted', label = f"Least Square Approx {ab_str}")

        ax.set_ylabel(f"Effective {ab_str}  [m^-1]")
        ax.set_xlabel("Wavelength [nm]")
    
        ax.grid()
        ax.legend()

        ## [Setting som etext with the reslting coefficients.]
#        ylim = ax.get_ylim()
#        xlim = ax.get_xlim()
#        ax.text(xlim*0.1, ylim*0.1, txt, bbox=props)
        ax.set_title(f'Constrained Least Squares Phytoplankton \n Community Estimation at Profile Date: {prof_sd}, Time: {prof_t}')
    
        fig.show()
    
    return 


def Plot_Spec_Lstsq_Hoffmuller_Time(N_profs, depthz, flort_profs, optaa_profs, phy_species, x_profs, phycolors, plot_residuals=False, residuals=None): 
    """
    This plots a Hoffmuller diagram of the ratios of different species opf phytoplankton in time

    It will be done by plotting a stacked bar plot.]
    """

    ## [The number of species.]
    N_phy = len(phy_species)

    ## [Init the bottom coordinate for the stacked bar plot.]
    bot = np.zeros(N_profs)

    ## [The colors for each species.]
#    colors = ['b','r','g','y','c','m','k','w'][:N_phy]
#    cmap = mpl.cm.get_cmap('Set1')

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
        ax_ls.bar(dt, x_profs[k,:], bottom=bot, color = phycolors[phy_species[k]], label=phy_species[k], width=1.0, alpha=.8)

        ## [Update the bottom of the bar plot.]
        bot = bot + x_profs[k,:]

    ax_ls.set_title(f'Constrained Least Squares \n Phytoplankton Community Estimation at Depth: {depthz}m')
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
        zi_flort = (abs(flort_profs[k]['depth'] - depthz)).argmin()
        chla[k] = flort_profs[k]['fluorometric_chlorophyll_a'].data[zi_flort]

    ## [Adding chla as twin plot to least squares.]
    ax_cl.plot(dt, chla, 'k', label='Chl-a', linewidth=2)
    ax_cl.set_ylim(bottom=0)
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
        

def Plot_Spec_Lstsq_Hoffmuller_Depth(N_depths, flort_prof, optaa_prof, phy_species, x_profs, phycolors, plot_residuals=False, residuals=None): 
    """
    This plots a Hoffmuller diagram of the ratios of different species opf phytoplankton at different depths for a single profile.

    It will be done by plotting a stacked bar plot.]
    """

    ## [The number of species.]
    N_phy = len(phy_species)

    ## [Init the bottom coordinate for the stacked bar plot.]
    bot = np.zeros(N_depths)

    ## [The colors for each species.]
#    colors = ['b','r','g','y','c','m','k','w'][:N_phy]
#    cmap = mpl.cm.get_cmap('Set1')

    ## [The true depth values.]
    depths = optaa_prof['depth'][:N_depths].data

    ## [The optaa depth coordinate]
    depthi = np.arange(N_depths)
    ## [If doing a bunch of depths, then just label every fifth depth.]
    if N_depths > 50:
        ytick_labels = []
        yticks = []
        for k in range(0,N_depths, int(N_depths/10)):
            ytick_labels.append(depths[k])
            yticks.append(k)
    else: 
        ytick_labels = depths
        yticks = depthi

    ## [Plot the residual on the same figure.]
    if plot_residuals: 
        fig, axs = plt.subplots(nrows=1, ncols=2)
        ## [Least square axis.]
        ax_ls = axs[0]
        ## [Residual axis.]
        ax_rs = axs[1]
    else:
        fig, ax = plt.subplots()
        ## [Least square axis.]
        ax_ls = ax
 
    ## [Loop over the species.]
    for k in range(N_phy): 
        
        ## [plot the phy species layer of the bar plot.]
#        ax_ls.barh(depthi, x_profs[k,:], left=bot, color = cmap(k/(2*N_phy)), label=phy_species[k], height=1.0, alpha=0.8)
        ax_ls.barh(depthi, x_profs[k,:], left=bot, color = phycolors[phy_species[k]], label=phy_species[k], height=1.0, alpha=0.8)

        ## [Update the bottom of the bar plot.]
        bot = bot + x_profs[k,:]

    ## [The profile start date for the title.]
    prof_sd = str(optaa_prof['time'].data.astype('datetime64[D]')[0])
    prof_t = str(optaa_prof['time'].data.astype('datetime64[m]')[0])[-5:]
    ax_ls.set_title(f'Constrained Least Squares Phytoplankton \n Community Estimation at Profile Date: {prof_sd}, Time: {prof_t}')
    ax_ls.set_xlabel('Fractional Concentration')
    ax_ls.set_ylabel('Depth [m]')
    ax_ls.set_yticks(yticks)
    ax_ls.set_yticklabels(ytick_labels)
    ax_ls.legend(loc=2)
    ax_ls.grid(axis='x')

    ## [Getting the chla at the givenz level in each profile.]
    chla = np.zeros(N_depths)
    for k in range(N_depths):
        ## [This aligns the chla concentration at a near z-coordinate of OPTAA]
        zi_flort = (abs(flort_prof['depth'] - depths[k])).argmin()
        chla[k] = flort_prof['fluorometric_chlorophyll_a'].data[zi_flort]

    ## [Chla twin axis on least squares plot.]
    ax_cl = ax_ls.twiny()
    ## [Adding chla as twin plot to least squares.]
    ax_cl.plot(chla, depthi, 'k', label='Chl-a', linewidth=2)
    ax_cl.set_xlim(left=0)
    ax_cl.set_xticks(np.linspace(ax_cl.get_xticks()[0], ax_cl.get_xticks()[-1], len(ax_ls.get_xticks())))
    ax_cl.set_xlabel(f"Fluorometirc Chl-a [{flort_prof['fluorometric_chlorophyll_a'].attrs['units'][0]}]")
    ax_cl.legend(loc=1)

    if plot_residuals: 
        ## [Plot the residuals.]
        ax_rs.barh(depthi, residuals, color = 'k', height=1.0)
        ax_rs.set_title('Two Norm of Residual Vector')
        ax_rs.set_xlabel(r'$\frac{|\mathbf{A} \mathbf{x} - \mathbf{y} |_2}{|\mathbf{y}|_2}$')
        ax_rs.set_yticks(yticks)
        ax_rs.set_yticklabels([])
        ax_rs.grid(axis='x')

    fig.show()

    return 
     
