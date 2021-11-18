"""
Created November 10 2021 

Author Miles MIller 
"""

import matplotlib.pyplot as plt
import numpy as np


def Plot_Comparison(ax, x, y, label, xlim=None, ylim=None, bf_line=False): 
    """
    Plots the given values on a given axes
    """

    ax.plot(x, y,'o', fillstyle='none', label=label, markersize=5)
    ax.plot(x, x, 'k')
    if bf_line: 
        m, b = np.polyfit(x,y,1) 
        ax.plot(x, m*x +b, 'r:')
    if xlim == None: 
        xlim = x.max()
    if ylim == None: 
        ylim = y.max()
    ax.set_xlim([0, xlim])
    ax.set_ylim([0, ylim])
    if bf_line: 
        return m, b 
    else:
        return ax 


def Plot_Frequency(ax, y, N_bins, label): 
    """ 
    Plots the frequency of the data as a line
    """

    hist, bin_edges = np.histogram(y, bins = N_bins)
   
    bin_centers = bin_edges[:-1] + (abs(bin_edges[0] - bin_edges[1]))

    ax.plot(bin_centers, hist, label)
  
    return ax
