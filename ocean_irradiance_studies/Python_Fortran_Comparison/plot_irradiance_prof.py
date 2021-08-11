# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 09:23:45 2021

@author: miles

ABOUT 
------

This is a file for looking at irradiance profile from NetCDF ROMS out. 
"""
## USER MOD
from ocean_irradiance_module.Read_ROMS_Out import ROMS_netcdf 
## GLOBAL MOD 
import numpy 
import matplotlib.pyplot as plt








if __name__ == '__main__':
    file = 'C:/Users/miles/RESEARCH/Current_Work/Ocean_Irradiance/ocean_irradiance_studies/Python_Fortran_Comparison/roms_his_phy.nc'
    R_nc = ROMS_netcdf(file,Init_ROMS_Irr_Params=(True))
    Ed0 = R_nc.Ed0

    
    ## Location and time step
    i = 10
    j = 10
    dt = 1
    ## Worth noting that the netcdf is read by python opposite of ROMS indexing
    z1 =  R_nc.z_irr1[dt,:,j,i]
    # Ed1 = R_nc.Ed1[dt,:,:,:]
    
    
    ## Plotting 
    # fig, ax = plt.subplots()
    # ax.plot(R_nc.Ed1[dt,:,j,i], z1, label = 'Ed1')
    # ax.plot(R_nc.Es1[dt,:,j,i], z1, label = 'Es1')
    # ax.plot(R_nc.Eu1[dt,:,j,i], z1, label = 'Eu1')
    
    # ax.set_title('Irradiance Profiles')
    # ax.legend()
    # ax.grid()
    # ax.set_ylabel('z [m]')
    # ax.set_xlabel('Irradiance')
    

