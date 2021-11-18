# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 10:24:36 2021

@author: miles

Interpolation Function 

    This scheme is to be tested on ROMS output and compared top the numpy module
    and then used in FORTRAN irradiance code. 
"""

import numpy as np
import os
cwd = os.getcwd()
os.chdir('../ocean_irradiance_module')
import Ocean_Irradiance as OI
os.chdir(cwd)
import matplotlib.pyplot as plt


def Interp(x,xd,yd):
    """
    Linear spline interpolation scheme. 

    Parameters
    ----------
    x : 1-D Array [A]
        The x-coordinates at which to evaluate the interpolated values. 
    xd : 1-D Array [B]
        The x-coordinates of the data points.
    yd : 1-D Array [B]
        The y-coordinates of the data points. 

    Returns
    -------
    y : 1-D Array [A]
        The interpolated values. 

    """
    
    N = len(xd)
    Nm1 = N-1
    y = np.zeros_like(x)
    p = 0 
    for k in range(Nm1):
        while xd[k] <= x[p] < xd[k+1]:
            y[p] = yd[k] + (x[p] - xd[k])*((yd[k+1] - yd[k])/(xd[k+1]-xd[k]))
            if p == len(x)-1:
                break
            else:
                p = p + 1
            
            
    return y

if __name__ == '__main__':
    zbot = -80

    xd = np.linspace(zbot,0,42)
    yd = xd**2
    
    zbot2 = -90
    x = OI.Log_Trans(zbot2,30)
    y = Interp(x, xd, yd)
    
    fig,ax = plt.subplots()
    ax.plot(xd,yd, 'ro', label='Data')
    ax.plot(x,y,'bo', label='Interp')
    ax.grid()
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
    