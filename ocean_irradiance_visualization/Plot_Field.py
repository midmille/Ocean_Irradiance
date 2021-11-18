"""
Created on Tue October 12 10:34:02 2021

@author: Miles Miller
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def Plot_Field(field, lat, lon, title, cbar_label, vmin=None, vmax=None, fig=None):

    cbar_shrink = 1

    fig = plt.figure()
    
    ## The first value plot
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
    #ax.gridlines()
    im = ax.pcolormesh(lon, lat, field, cmap='nipy_spectral', 
                     transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin )
    fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = cbar_label)
    ax.set_title(title)  
    ylims = ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
    ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))

    fig.show()
 
    return 


def Plot_Fields(subplot_rows, subplot_cols, subplot_loc, field, lat, lon, title, cbar_label, fig=None, vmin=None, vmax=None):

    cbar_shrink = 1

    if fig == None:
        fig = plt.figure() 

    ## The first value plot
    ax = fig.add_subplot(subplot_rows, subplot_cols, subplot_loc, projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
    #ax.gridlines()
    im = ax.pcolormesh(lon, lat, field, cmap='nipy_spectral', 
                     transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin )
    fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = cbar_label)
    ax.set_title(title)  
    ylims = ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
    ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))
 
    return fig


def Plot_Scatter(field, lat, lon, title, cbar_label, vmin=None, vmax=None, fig=None):

    cbar_shrink = 1

    fig = plt.figure()
    
    ## The first value plot
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND, color='grey', alpha=.5)
    #ax.gridlines()
    s = 15
    im = ax.scatter(lon, lat, c=field, s=s, cmap='nipy_spectral', 
                     transform=ccrs.PlateCarree(), vmax=vmax, vmin=vmin )
    fig.colorbar(im, ax=ax, shrink=cbar_shrink, label = cbar_label)
    ax.set_title(title)  
    ylims = ax.set_ylim(ymin=np.min(lat), ymax=np.max(lat))
    ax.set_xlim(xmin=np.min(lon), xmax=np.max(lon))

    fig.show()
 
    return 







#if __name__ == '__main__': 

#    import argparse 

#    parser = argparse.ArgumentParser(description='Ocean Irradiance Plotting Toolkit')
#    parser.add_argument('field', help = "Field to be plotted" )
#    parser.add_argument('lat', help = "Lattitude Array" )
#    parser.add_argument('lon', help = "Longtitude Array" )
#    parser.add_argument('title', help = "" )
#    parser.add_argument('--plot', action='store_true', help="Visualization of Result")
#    args = parser.parse_args()
 
   

