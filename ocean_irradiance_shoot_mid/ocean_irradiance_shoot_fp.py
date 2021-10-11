# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 10:23:00 2021

@author: miles

This is the implementation of using the shoot method to a fitting point to solve 
the three stream irradiance bvp. 

The algorithim will follow closely that as described by Numerical Recipes
Third Edition chapter 18.2.
"""

import ocean_irradiance_module.Ocean_Irradiance as OI
import numpy as np



def Irradiance_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, r_s, r_u, v_d, v_s, v_u, 
                   direction = 'up'):
    
    
    """
    The RK4 algorithim that will compute iterate down from the surface of the water r
    column given an initial value. The initial value should be the very last index in the array.

    
    Assumes that ROMS orientation of veritcal grid. 
    
    keyword arg 'direction' implies the direction of the iteration:
        'down' ==> from surface to bottom
        'up' ==> from bottom to surface
    """
    
    if direction == 'down': 

        for k in range(Nm1-1 , -1, -1) :
    
            dz = z[k] - z[k+1]
            dzo2 = dz / 2
            dzo2o2 = dzo2 / 2;
              
            dEddz1 = c_d[k]*Ed[k+1]
            dEddz2 = c_d[k]*(Ed[k+1]+(dEddz1*dzo2))
            dEddz3 = c_d[k]*(Ed[k+1]+(dEddz2*dzo2))
            dEddz4 = c_d[k]*(Ed[k+1]+(dEddz3*dz))
            Ed[k] = Ed[k+1] + (( (dEddz1/6)+(dEddz2/3)+(dEddz3/3)+(dEddz4/6) )*dz)
              
            ## to get Edmid, integrate Ed eq only down to mid-point
            dEddz2 = c_d[k]*(Ed[k+1]+(dEddz1*dzo2o2))
            dEddz3 = c_d[k]*(Ed[k+1]+(dEddz2*dzo2o2))
            dEddz4 = c_d[k]*(Ed[k+1]+(dEddz3*dzo2))
            Edmid = Ed[k+1] + (( (dEddz1/6)+(dEddz2/3)+(dEddz3/3)+(dEddz4/6) )*dzo2)
              
            dEsdz1 = - (b_f[k]/v_d)*Edmid + ((a[k]+r_s*b_b[k])/v_s)*Es[k+1] - (r_u*b_b[k]/v_u)*Eu[k+1]
               
            dEudz1 =  (b_b[k]/v_d)*Edmid + (r_s*b_b[k]/v_s)*Es[k+1] - ((a[k]+r_u*b_b[k])/v_u)*Eu[k+1]
              
            dEsdz2 = - (b_f[k]/v_d)*Edmid + ((a[k]+r_s*b_b[k])/v_s)*((Es[k+1])+(dEsdz1*dzo2)) \
            - (r_u*b_b[k]/v_u)*(Eu[k+1]+(dEudz1*dzo2))
              
            dEudz2 =  (b_b[k]/v_d)*Edmid + (r_s*b_b[k]/v_s)*(Es[k+1]+(dEsdz1*dzo2)) \
            - ((a[k]+r_u*b_b[k])/v_u)*(Eu[k+1]+(dEudz1*dzo2))
              
            dEsdz3 = - (b_f[k]/v_d)*Edmid + ((a[k]+r_s*b_b[k])/v_s)*(Es[k+1]+(dEsdz2*dzo2)) \
            - (r_u*b_b[k]/v_u)*(Eu[k+1]+(dEudz2*dzo2))
              
            dEudz3 =  (b_b[k]/v_d)*Edmid + (r_s*b_b[k]/v_s)*(Es[k+1]+(dEsdz2*dzo2)) \
            - ((a[k]+r_u*b_b[k])/v_u)*(Eu[k+1]+(dEudz2*dzo2))
              
            dEsdz4 = - (b_f[k]/v_d)*Edmid + ((a[k]+r_s*b_b[k])/v_s)*(Es[k+1]+(dEsdz3*dz)) \
            - (r_u*b_b[k]/v_u)*(Eu[k+1]+(dEudz3*dz))
              
            dEudz4 =  (b_b[k]/v_d)*Edmid + (r_s*b_b[k]/v_s)*(Es[k+1]+(dEsdz3*dz)) \
            - ((a[k]+r_u*b_b[k])/v_u)*(Eu[k+1]+(dEudz3*dz))
              
              ## RK4
            Es[k] = Es[k+1] + (( (dEsdz1/6)+(dEsdz2/3)+(dEsdz3/3)+(dEsdz4/6))*dz)
            Eu[k] = Eu[k+1] + (( (dEudz1/6)+(dEudz2/3)+(dEudz3/3)+(dEudz4/6))*dz)
    
    elif direction == 'up':

        Edmid = np.zeros(Nm1) ##cell centers
        ## integrate Es down the water column
        for k in range(Nm1-1, -1, -1) :
         # carry out a RK4 algorithm
    
         # Ackelson-Dutkiewicz equations
         # z is +ve upward, opposite to that in Dutkiewicz
         # steping DOWNWARD through water column means dz should be less than 0
            dz = z[k] - z[k+1]
            dzo2 = dz / 2
            dzo2o2 = dzo2 / 2;
    
         ## to get Edmid, integrate Ed eq only down to mid-point
            dEddz1 = c_d[k]*Ed[k+1]
            dEddz2 = c_d[k]*(Ed[k+1]+(dEddz1*dzo2o2))
            dEddz3 = c_d[k]*(Ed[k+1]+(dEddz2*dzo2o2))
            dEddz4 = c_d[k]*(Ed[k+1]+(dEddz3*dzo2))
            Edmid[k] = Ed[k+1] + (( (dEddz1/6)+(dEddz2/3)+(dEddz3/3)+(dEddz4/6) )*dzo2)
            
            dEddz1 = c_d[k]*Edmid[k]
            dEddz2 = c_d[k]*(Edmid[k]+(dEddz1*dzo2o2))
            dEddz3 = c_d[k]*(Edmid[k]+(dEddz2*dzo2o2))
            dEddz4 = c_d[k]*(Edmid[k]+(dEddz3*dzo2))
            Ed[k] = Edmid[k] + (( (dEddz1/6)+(dEddz2/3)+(dEddz3/3)+(dEddz4/6) )*dzo2)

        for k in range(Nm1) :
            

            dz = - z[k] + z[k+1]
            dzo2 = dz / 2
            dzo2o2 = dzo2 / 2;
            

            dEsdz1 = - (b_f[k]/v_d)*Edmid[k] + ((a[k]+r_s*b_b[k])/v_s)*Es[k] - (r_u*b_b[k]/v_u)*Eu[k]
         
            dEudz1 =  (b_b[k]/v_d)*Edmid[k] + (r_s*b_b[k]/v_s)*Es[k] - ((a[k]+r_u*b_b[k])/v_u)*Eu[k]
    
            dEsdz2 = - (b_f[k]/v_d)*Edmid[k] + ((a[k]+r_s*b_b[k])/v_s)*((Es[k])+(dEsdz1*dzo2)) \
            - (r_u*b_b[k]/v_u)*(Eu[k]+(dEudz1*dzo2))
            
            dEudz2 =  (b_b[k]/v_d)*Edmid[k] + (r_s*b_b[k]/v_s)*(Es[k]+(dEsdz1*dzo2)) \
            - ((a[k]+r_u*b_b[k])/v_u)*(Eu[k]+(dEudz1*dzo2))
    
            dEsdz3 = - (b_f[k]/v_d)*Edmid[k] + ((a[k]+r_s*b_b[k])/v_s)*(Es[k]+(dEsdz2*dzo2)) \
            - (r_u*b_b[k]/v_u)*(Eu[k]+(dEudz2*dzo2))
            
            dEudz3 =  (b_b[k]/v_d)*Edmid[k] + (r_s*b_b[k]/v_s)*(Es[k]+(dEsdz2*dzo2)) \
            - ((a[k]+r_u*b_b[k])/v_u)*(Eu[k]+(dEudz2*dzo2))
    
            dEsdz4 = - (b_f[k]/v_d)*Edmid[k] + ((a[k]+r_s*b_b[k])/v_s)*(Es[k]+(dEsdz3*dz)) \
            - (r_u*b_b[k]/v_u)*(Eu[k]+(dEudz3*dz))
            
            dEudz4 =  (b_b[k]/v_d)*Edmid[k] + (r_s*b_b[k]/v_s)*(Es[k]+(dEsdz3*dz)) \
            - ((a[k]+r_u*b_b[k])/v_u)*(Eu[k]+(dEudz3*dz))
    
            ## Euler
            # Es[k+1] = Es[k] + (dEsdz1*dz)
            # Eu[k+1] = Eu[k] + (dEudz1*dz)
    
            ## RK4
            Es[k+1] = Es[k] + (( (dEsdz1/6)+(dEsdz2/3)+(dEsdz3/3)+(dEsdz4/6))*dz)
            Eu[k+1] = Eu[k] + (( (dEudz1/6)+(dEudz2/3)+(dEudz3/3)+(dEudz4/6))*dz)
            

    return Ed, Es, Eu


def ocean_irradiance_shoot_up(hbot, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, N = 30, 
                     pt1_perc_zbot = True):
    
    
    """
    The main ocean_irradiance function that calculates the three stream model solution 
    following the equations and coefficients of Dutkiewicz (2015) and solved as a boundary 
    value problem using the shooting method. 

    Parameters
    ----------
    hbot : Float  
        True bottom depth of water column. 
    E_d_0 : Float
        Initial value of downward direct irradiance. 
    E_s_0 : Float
        Initial value of downward diffuse irradiance. 
    E_u_h : Float
        Boundary Condition on upwelling irradiance at h. 
    ab_wat : Tuple, length==2, (a_wat,b_wat).  
        Absorbtion and scattering coefficients for water. 
    coeffcients : Tuple, length == 5. 
        Coefficients taken from Dutkiewicz such as the average of cosines and sines. 
    phy : optional, default is None, else Phy object. 
        Gives information as according to Phy class on phytoplankton profile(s), 
        corresponding z-grid, and coefficients of absorbtion and scattering for each
        respective species of phytoplankton. 
    N : Float, default is 30
        The number of layers in the logarithmic grid. 
    pt1_perc_zbot : Boolean, default is True
        True refers to using the .1% light level as the zbot so long as that the magnitude 
        of the .1% light level is smaller than the magnitude of hbot. False refers
        to just using given hbot as zbot. 

    Returns
    -------
    Ed : 1-D Array
        Downward direct irradiance. 
    Es : 1-D Array
        Downward diffuse irradiance. 
    Eu : 1-D Array
        Downward diffuse irradiance. 
    z : 1-D Array
        The grid that the irradiances are calculated on. 

    """
    ## PARAMS FROM DUTKIEWICZ 2015 
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    ##N centers
    Nm1 = N - 1  
    
    ##initial guess... doesn't matter too much
    init_guess = .001
    # init_guess = 0
    
    Ed1 = np.full(N, init_guess)
    Es1 = np.full(N, init_guess)
    Eu1 = np.full(N, init_guess) 
    
    ## BCs included in the initial guess arrays
    Ed1[Nm1] = Ed0 ##surface 
    Es1[Nm1] = Es0 ##surface 
    Eu1[0] = Euh ##bottom
    
    ## Default number of shots for BVP shoot method solution
    shots =3
    
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

    
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        c_wat = (a_wat + b_wat)/v_d
        zbot_pt1perc = OI.zbot_func(Ed0, c_wat, z_phy)
        print(zbot_pt1perc)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    # z = Log_Trans(zbot, N) 
    ## linear z 
    z = np.linspace(zbot, 0, N)
    
    
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    if phy: 
        a = np.interp(z,z_phy,a)
        b = np.interp(z,z_phy,b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        
    ##coefficient of downward direct irradiance 
    c_d = (a+b)/v_d 
    # if N != len(Ed1) or N !=len(a)+1 or N !=len(b)+1 :
    #     print('lengths of z and Ed must be the same, and a&b should be 1 less.')
        

    b_b = .551*b
    b_f = b - b_b 
    

        
    Ed=np.copy(Ed1)
    Es=np.copy(Es1)
    Eu=np.copy(Eu1)
 
    Es0_tried = []
    Fmetric = []
     
    for jsh in range(shots) :
    # Integrate down from the top to ensure Ed(1) and Es(1) are good.
         if jsh == 0:
             dEs = 0 #-Eu[Nm1]
         elif jsh == 1:
         # for the first case, need some adjustment to get gradient.
             # dEs = max(0.01,0.03*Es[Nm1])
             dEs = 0.01
         else: 

             Jslope = (Fmetric[jsh-2]-Fmetric[jsh-1]) / (Es0_tried[jsh-2]-Es0_tried[jsh-1]) 
             dEs = -Fmetric[jsh-1]/Jslope
             # print(Jslope)
         # dEs = max(-Es[0], dEs)  # make sure Es can't go below 0
         # dEs = min(1-Es[0],dEs) # make sure Es can't go above 1
         # print(dEs)
         Es[0] = Es[0] + dEs
         Es0_tried.append(Es[0])
 
         # Edmid = np.zeros(Nm1)
         ## integrate Es down the water column
         ## The range does not actually go to k=-1, only to k=0. 
         ## i.e. does not include stop point of range. 
         # for k in range(Nm1-1 , -1, -1) :
         Ed, Es, Eu = Irradiance_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, 
                                     r_s, r_u, v_d, v_s, v_u, direction='up')

         Fmetric.append(Es0 - Es[Nm1])
     # Eu[-1] = Eu[-1] - .0001
     # Ed, Es, Eu = Irradiance_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, 
     #                             r_s, r_u, v_d, v_s, v_u)


    return Ed, Es, Eu, z



def ocean_irradiance_shoot_fp(hbot, fp, fpi, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, N = 30, 
                     pt1_perc_zbot = True):
    
    
    """
    The main ocean_irradiance function that calculates the three stream model solution 
    following the equations and coefficients of Dutkiewicz (2015) and solved as a boundary 
    value problem using the shooting method. This version shoots to a fixed point from both directions
    and solves for continuity. 

    Parameters
    ----------
    hbot : Float  
        True bottom depth of water column. 
    E_d_0 : Float
        Initial value of downward direct irradiance. 
    E_s_0 : Float
        Initial value of downward diffuse irradiance. 
    E_u_h : Float
        Boundary Condition on upwelling irradiance at h. 
    ab_wat : Tuple, length==2, (a_wat,b_wat).  
        Absorbtion and scattering coefficients for water. 
    coeffcients : Tuple, length == 5. 
        Coefficients taken from Dutkiewicz such as the average of cosines and sines. 
    phy : optional, default is None, else Phy object. 
        Gives information as according to Phy class on phytoplankton profile(s), 
        corresponding z-grid, and coefficients of absorbtion and scattering for each
        respective species of phytoplankton. 
    N : Float, default is 30
        The number of layers in the logarithmic grid. 
    pt1_perc_zbot : Boolean, default is True
        True refers to using the .1% light level as the zbot so long as that the magnitude 
        of the .1% light level is smaller than the magnitude of hbot. False refers
        to just using given hbot as zbot. 

    Returns
    -------
    Ed : 1-D Array
        Downward direct irradiance. 
    Es : 1-D Array
        Downward diffuse irradiance. 
    Eu : 1-D Array
        Downward diffuse irradiance. 
    z : 1-D Array
        The grid that the irradiances are calculated on. 

    """
    ## PARAMS FROM DUTKIEWICZ 2015 
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    ##N centers
    Nm1 = N - 1  
    
    ##initial guess... doesn't matter too much
    init_guess = .001
    # init_guess = 0
    
    Ed1 = np.full(N, init_guess)
    Es1 = np.full(N, init_guess)
    Eu1 = np.full(N, init_guess) 
    
    ## BCs included in the initial guess arrays
    Ed1[Nm1] = Ed0 ##surface 
    Es1[Nm1] = Es0 ##surface 
    Eu1[0] = Euh ##bottom
    
    ## Default number of shots for BVP shoot method solution
    shots = 20
    
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

    
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        c_wat = (a_wat + b_wat)/v_d
        zbot_pt1perc = OI.zbot_func(Ed0, c_wat, z_phy)
        print(zbot_pt1perc)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    # z = Log_Trans(zbot, N) 
    ## linear z 
    z = np.linspace(zbot, 0, N)
    ## The fpi is the 80% point, closer to surface
    fpi = int(.99 * N)
    fp = z[fpi]
    
    
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    if phy: 
        a = np.interp(z,z_phy,a)
        b = np.interp(z,z_phy,b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        
    ##coefficient of downward direct irradiance 
    c_d = (a+b)/v_d 
    # if N != len(Ed1) or N !=len(a)+1 or N !=len(b)+1 :
    #     print('lengths of z and Ed must be the same, and a&b should be 1 less.')
        

    b_b = .551*b
    b_f = b - b_b 

    
    ## Irradiances Up
    Edd=np.copy(Ed1)
    Esd=np.copy(Es1)
    Eud=np.copy(Eu1)
    
    ## Irradiances down
    Edu=np.copy(Ed1)
    Esu=np.copy(Es1)
    Euu=np.copy(Eu1)
 
    Es0_tried = []
    Eu0_tried = []
    
    Fmetric_Es = []
    Fmetric_Eu = []
    Fmetric = []
     
    for jsh in range(shots) :
    # Integrate down from the top to ensure Ed(1) and Es(1) are good.
         if jsh == 0:
             dEs = 0 #-Eu[Nm1]
             dEu = 0
         elif jsh == 1:
         # for the first case, need some adjustment to get gradient.
             dEu = .02 #max(0.01,0.03*Es[Nm1])
             dEs = 0.01
         else: 

             Jslope_Es = (Fmetric_Es[jsh-2]-Fmetric_Es[jsh-1]) / (Es0_tried[jsh-2]-Es0_tried[jsh-1]) 
             Jslope_Eu = (Fmetric_Eu[jsh-2]-Fmetric_Eu[jsh-1]) / (Eu0_tried[jsh-2]-Eu0_tried[jsh-1])
             
             dEs = -Fmetric_Es[jsh-1]/Jslope_Es
             dEu = -Fmetric_Eu[jsh-1]/Jslope_Eu
             # print(Jslope)
         # dEs = max(-Es[0], dEs)  # make sure Es can't go below 0
         # dEs = min(1-Es[0],dEs) # make sure Es can't go above 1
         # print(dEs)
         Esu[0] = Esu[0] + dEs
         Es0_tried.append(Esu[0])
         
         Eud[Nm1] = Eud[Nm1] + dEu
         Eu0_tried.append(Eud[Nm1])
 
         # Edmid = np.zeros(Nm1)
         ## integrate Es down the water column
         ## The range does not actually go to k=-1, only to k=0. 
         ## i.e. does not include stop point of range. 
         # for k in range(Nm1-1 , -1, -1) :
         Edu, Esu, Euu = Irradiance_RK4(Nm1, Edu, Esu, Euu, z, a, b, c_d, b_b, b_f, 
                                     r_s, r_u, v_d, v_s, v_u, direction='up')
         Edd, Esd, Eud = Irradiance_RK4(Nm1, Edd, Esd, Eud, z, a, b, c_d, b_b, b_f, 
                                     r_s, r_u, v_d, v_s, v_u, direction='down')

         # Fmetric.append(Es0 - Es[Nm1])
         Fmetric_Es.append(Esu[fpi] - Esd[fpi])
         Fmetric_Eu.append(Euu[fpi] - Eud[fpi])
         # Fmetric.append(np.sqrt((Esu[fpi] - Esd[fpi])**2 + (Euu[fpi] - Eud[fpi])**2))

     # Eu[-1] = Eu[-1] - .0001
     # Ed, Es, Eu = Irradiance_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, 
     #                             r_s, r_u, v_d, v_s, v_u)

    Ed = np.append(Edu[:fpi], Edd[fpi:])
    Es = np.append(Esu[:fpi], Esd[fpi:])
    Eu = np.append(Euu[:fpi], Eud[fpi:])
    return Ed, Es, Eu, z, fpi


def Demo(method='shoot_up'): 
    ## Demonstration of finding the irradiance profiles for an artificial 
    ## phytoplankton concentration profile
    
    ## Absorbtion and scattering coefficients are taken from Dutkiwicz 2015 
    ## for the Diatom species. 
    import matplotlib.pyplot as plt
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    PI = Param_Init()
    
    N = 300
    Nm1 = N-1 
    lam =443
    
    z = np.linspace(-600,0,N)

    phy_prof = OI.artificial_phy_prof(z, -5, 10,.5)
    # ROMS_point = np.genfromtxt('ChrisData_good_point.csv', delimiter=',')
    # phy_prof = ROMS_point[1:,2]
    # print(phy_prof)
    # z = ROMS_point[1:,0]
    # phy_prof = np.full(len(z), 1)
    
    ab_wat = abscat(lam, 'water')
    
    a_phy, b_phy = abscat(lam, 'Syn')
    
    
    ## Define the Phytoplankton class.
    phy = OI.Phy(z, phy_prof, a_phy, b_phy)


    ## The fixed point position: 
    fp = -50
    fpi =0

    zbot = z[0]
    
    if method == 'shoot_up':
        
        Ed, Es, Eu, zarr = ocean_irradiance_shoot_up(zbot,PI.Ed0,PI.Es0,PI.Euh,ab_wat, PI.coefficients, 
                                            phy=phy, N=N, pt1_perc_zbot = True)
     
    if method == 'shoot_down':
        Ed, Es, Eu, zarr = OI.ocean_irradiance(zbot, PI.Ed0, PI.Es0, PI.Euh, ab_wat,  PI.coefficients,
                                               phy=phy, N=N, pt1_perc_zbot=True)
        
    if method == 'shoot_fp': 
        Ed, Es, Eu, zarr, fpi = ocean_irradiance_shoot_fp(zbot, fp, fpi, PI.Ed0, PI.Es0, PI.Euh, 
                                                     ab_wat, PI.coefficients, phy = phy, N = N, 
                                                     pt1_perc_zbot = True)
        
    if method == 'scipy':
        Ed, Es, Eu, zarr = OI.ocean_irradiance(zbot, PI.Ed0, PI.Es0, PI.Euh, ab_wat,  PI.coefficients,
                                               phy=phy, N=N, pt1_perc_zbot=False, use_bvp_solver=True)    

    if method == 'dut':
        Ed, Es, Eu, zarr = OI.ocean_irradiance_dutkiewicz(zbot, PI.Ed0, PI.Es0, PI.Euh, ab_wat,  
                                                          PI.coefficients, phy=phy, N=N, pt1_perc_zbot=False)
    ## Plotting the Results
    #-------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    
    ax1.plot(phy_prof, z)
    if method == 'shoot_fp':
        ax1.hlines( zarr[fpi], min(phy_prof),max(phy_prof), color='r')
    ax1.set_xlabel('Concentration')
    ax1.set_ylabel('Z [m]')
    ax1.set_title('Phytoplankton Concentration Profile')
    ax1.grid()
    
    ax2.plot(Ed, zarr, label='Ed')
    ax2.plot(Es, zarr, label='Es')
    ax2.plot(Eu, zarr, label='Eu')
    ax2.set_xlim(-0.1,1)
    ax2.set_xlabel('Irradiance')
    ax2.set_title('Resulting Irradiance Profiles')
    ax2.legend()
    ax2.grid()
    
    plt.show()
    return zarr, Ed, fpi




if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Ocean Irradiance Fucntion')

    parser.add_argument('--demo', action='store_true', 
                        help='Run the demo and plot')
    
    args = parser.parse_args()
    
    # if args.demo: 
    zarr, Ed, fpi = Demo('shoot_fp')
