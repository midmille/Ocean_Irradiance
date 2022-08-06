"""
Created November 10 2021 

Author Miles MIller 
"""

import matplotlib.pyplot as plt
import numpy as np


def Plot_Comparison(ax, x, y, title, label, xlabel, ylabel, xlim=None, ylim=None, marker='o', color=None, title_fontsize=10, markersize=5, linewidth=1): 
    """
    Plots the given values on a given axes
    """
    ## [remove any 'nans': 
    nonan = ~np.isnan(x) * ~np.isnan(y)
    x = x[nonan]
    y = y[nonan]

    if color==None:
        ax.plot(x, y,'o', fillstyle='none', label=label, markersize=markersize, linewidth=linewidth)
    else:
        ax.plot(x, y,'o', fillstyle='none', markeredgecolor=color, label=label, markersize=markersize, linewidth=linewidth)
    if xlim == None: 
        ax.relim()
    elif ylim == None: 
        ax.relim()
    else:
        ax.set_xlim([0, xlim])
        ax.set_ylim([0, ylim])
    ## Using the limits as the coordinates for the x=y line.
    ylims = ax.get_ylim()
    xlims = ax.get_xlim()
    low_lim_min = min(ylims[0], xlims[0]) 
    high_lim_max = max(ylims[1], xlims[1])
    ax.plot([low_lim_min, high_lim_max], [low_lim_min, high_lim_max], 'k')

    ax.set_title(title, fontsize =title_fontsize)
    ax.set_xlabel(xlabel, fontsize =title_fontsize)
    ax.set_ylabel(ylabel, fontsize =title_fontsize)

    ## [Calculating and printing correlation statisitcs.]
    RMS = np.sqrt(np.mean(((y-x)/x)**2))
    mean_bias = np.mean(y-x)
    mean_ratio = np.mean(y/x)
    slope, intercept = np.polyfit(x,y,1)
    N = len(x)

    print("RMSRE:", RMS)
    print("mean bias:", mean_bias)
    print("mean_ratio:", mean_ratio)
    print("slope:", slope)
    print("intercept:", intercept)
    print("N:", N)
    print(f'{round(RMS,4)} & {round(mean_bias,4)} &{round(mean_ratio,4)}&{round(slope,4)}& {round(intercept,4)}& {N}')
    
    return ax 


#def Plot_Comparison(ax, x, y, label, xlim=None, ylim=None, bf_line=False): 
    #"""
    #Plots the given values on a given axes
    #"""
#
    #ax.plot(x, y,'o', fillstyle='none', label=label, markersize=5)
    #ax.plot(x, x, 'k')
    #if bf_line: 
        #m, b = np.polyfit(x,y,1) 
        #ax.plot(x, m*x +b, 'r:')
    #if xlim == None: 
        #xlim = x.max()
    #if ylim == None: 
        #ylim = y.max()
    #ax.set_xlim([0, xlim])
    #ax.set_ylim([0, ylim])
    #if bf_line: 
        #return m, b 
    #else:
        #return ax 


def Plot_Frequency(ax, y, N_bins, label, bin_edges=[]): 
    """ 
    Plots the frequency of the data as a line
    """


    y = y[~np.isnan(y)]
    if len(bin_edges) > 0:
        hist, bin_edges = np.histogram(y, bins = bin_edges)
    else:
        hist, bin_edges = np.histogram(y, bins = N_bins)
   
    bin_centers = bin_edges[:-1] + (abs(bin_edges[0] - bin_edges[1]))

    ax.plot(bin_centers, hist, label =label)
  
    return ax, bin_edges
