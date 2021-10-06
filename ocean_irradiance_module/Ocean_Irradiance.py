"""
Created on Thu Dec 10 08:52:28 2020

@author: Miles Miller
"""
import numpy as np
import scipy as sp 
from scipy import integrate



class Phy:
    """
    The Phytoplankton Class includes:
        --> z coordinates corresponding to phytoplankton concentrations [m (< 0)]
                -- z should be a column vector
        --> Phytoplankton chlorophyll concentrations [mg chl-a m^-3]
                --2-D array with columns corresponding to diff. phytoplankton
        --> absorbtion coefficient, a [m^-1]
                -- Vector corresponding to different phytoplankton
        --> backscatter coefficient, b [m^-1]
                -- Vector corresponding to different phytoplankton
    """
    def __init__(self, z, phy, a, b):
        self.z = z
        self.phy = phy
        self.a = a
        self.b = b
        self.Nz = len(z)
        if phy.ndim == 1:
            self.Nphy = 1
        elif phy.ndim == 2:
            self.Nphy = np.shape(phy)[1]
            assert self.Nphy == len(a) and self.Nphy == len(b)

        assert np.shape(phy)[0] == self.Nz
        
        return 


def analytical_Ed(zarr,c, Ed0): 
    
    
    """
    Parameters
    ----------
    zarr : vertical 1-D array
        The vertical array of the depth from 0 --> negative depth.

    Returns
    -------
    Downward Direct Irradiance
    

    """
   
    Ed = Ed0*np.exp(c*zarr)
    return Ed 


def numerical_Ed(z, c_d, Ed0):
    
    N = len(z)
    Nm1 = N - 1
    
    Ed = np.zeros(N)
    Ed[Nm1] = Ed0

    for k in range(Nm1-1, -1, -1) :
        
        dz = z[k] - z[k+1]
        dzo2 = dz / 2
        
        dEddz1 = c_d[k]*Ed[k+1]
        dEddz2 = c_d[k]*(Ed[k+1]+(dEddz1*dzo2))
        dEddz3 = c_d[k]*(Ed[k+1]+(dEddz2*dzo2))
        dEddz4 = c_d[k]*(Ed[k+1]+(dEddz3*dz))
        Ed[k] = Ed[k+1] + (( (dEddz1/6)+(dEddz2/3)+(dEddz3/3)+(dEddz4/6) )*dz)
        
    return Ed


def numerical_Ed_2(z, c, Ed0):
    
    N = len(z)
    Nm1 = N-1
    Nm2 = N-2
    #assert N == len(c)
    A = np.zeros((N,N))
    dz = abs(zarr[0] - zarr[1])
    m = 1/(2*dz)
    row_index = [i for i in range(Nm2,0,-1)]
    j = -1
    for i in row_index: 
        A[i,j] = m
        A[i,j-1] = -c[i]
        A[i,j-2] = -m
        j -= 1 
    A[-1, -1] = 1
    A[0,0] = -(c[0]+2*m)
    A[0,1] = 2*m 
    B = np.zeros(N) #must have same number of rows as A 
    B[-1] = Ed0 ##setting the first number to the amplitude at surface
    lu, piv = sp.linalg.lu_factor(A)
    Ed = sp.linalg.lu_solve((lu, piv), B)
    """
    At = np.transpose(A)
    AtA = At @ A
    AtA_inverse = np.linalg.inv(AtA)
    AtB = At @ B
    x = (AtA_inverse @ AtB)
    """
    return Ed


def zbot_func(E_d_0, c, z):
    """
    Finds the zbot for at which light ha attenuated to .1% of its surface value 
    for water only coeffients

    Parameters
    ----------
    E_d_0 : Float 
        Initial Value for E_d. 
    a_wat : Float 
        Absorbtion coefficient of water. 
    b_wat : Float
        Scattering Coefficient of water. 

    Returns
    -------
    zbot : Float 
        .01% light level zbot. 

    """
    zbots = np.linspace(0, -1000, 2001) 
    # c = (a_wat + b_wat) / v_d
#    c = np.flip(c)
#    z = np.flip(z)
    Ed = analytical_Ed(zbots, c, E_d_0)
    for k, Ed_i in enumerate(Ed) :
        EdoE_d_0 = Ed_i / E_d_0
        if EdoE_d_0 < .01 :
            zbot = zbots[k] 
            return zbot
   
        
def Log_Trans(zbot,Nlayers):
    """
    Constructs a Log Transformed z grid. 

    Parameters
    ----------
    zbot : Float
        zbot at .01% light level.
    Nlayers : Int 
        Number of layers to create log grid. 

    Returns
    -------
    zarr : 1-D Array 
        The resulting z-grid. 

    """
    zarr_exp = np.linspace(np.log(abs(zbot)),0, Nlayers)
    zarr = -np.exp(zarr_exp)
    
    return zarr


def ocean_irradiance_scipy(zarr, Ed0, Es0, Euh, a, b, coefficients): 
    """
    This fucntion uses the scipy.integrate.solve_bvp function to solve the ocean 
    irradiance boundary value problem. This is normally taken to be the truth in 
    sensitivity studies. 

    Parameters
    ----------
    zarr : 1-D Array [N]
        The vertical z coordinates. 
    Ed0 : Scalar
        Surface boundary value of Ed. 
    Es0 : Scalar 
        Surface boundary value of Es. 
    Euh : Scalar 
        Bottom boundary value of Eu.
    a : 1-D Array [N]
        The total absorbtion coefficient array. 
    b : 1-D Array[N]
        The total scattering coefficient array. 
    coefficients : Tuple [5]
        The coefficients taken from Dutkiewicz 2015. Averaging of cosines and such.

    Returns
    -------
    Ed, Es, Eu : 1-D Arrays [N]
        Solution to ocean irradiance. 

    """

    
    def derivEdz(z,E):
        
        E_d = E[0,:]
        E_s = E[1,:]
        E_u = E[2,:]
        
        a_r = np.interp(z,zarr,a)
        b_r = np.interp(z,zarr,b)
        
        b_b = .551*b_r 
        b_f = b_r - b_b 
        
        ## PARAMS FROM DUTKIEWICZ 2015 
        r_s, r_u, v_d, v_s, v_u = coefficients
        
        dEddz = (a_r+b_r)/v_d*E_d
        dEsdz = -b_f/v_d*E_d +(a_r+r_s*b_b)/v_s*E_s    - r_u*b_b/v_u*E_u
        dEudz =  b_b/v_d*E_d    + r_s*b_b/v_s*E_s - (a_r+r_u*b_b)/v_u*E_u
        
        dEdz = np.array([dEddz, dEsdz, dEudz])
        
        return dEdz
        
    def Ebcs(E_at_h, E_at_0):

        
        return np.array([E_at_0[0] - Ed0, E_at_0[1] - Es0, E_at_h[2] - Euh])
    
    N = len(zarr)
    
    Eguess = np.full((3, N), 1.0)

    res = integrate.solve_bvp(derivEdz, Ebcs, zarr, Eguess)

    y = res.y

    Ed = y[0]
    Es = y[1]
    Eu = y[2]
    
    return Ed, Es, Eu 


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


def Scipy_RK4(Nm1, Ed, Es, Eu, zarr, a, b, c_d, b_b, b_f, r_s, r_u, v_d, v_s, v_u): 
    
    """
    This follows a similiar format of the Irradiance_RK4 routine but computes the 
    initial value problem using the scipy RK45 algorithim. 
    """
    
    
    def derivEdz(z, E):
    
        E_d, E_s, E_u = E
        
        a_r = np.interp(z,zarr,a)
        b_r = np.interp(z,zarr,b)
        
        
        b_b = .551*b_r 
        b_f = b_r - b_b 
        
        dEddz = (a_r+b_r)/v_d*E_d
        dEsdz = -b_f/v_d*E_d +(a_r+r_s*b_b)/v_s*E_s    - r_u*b_b/v_u*E_u
        dEudz =  b_b/v_d*E_d    + r_s*b_b/v_s*E_s - (a_r+r_u*b_b)/v_u*E_u
        
        dEdz = np.array([dEddz, dEsdz, dEudz])
        
        return dEdz
    
    res = integrate.solve_ivp(derivEdz, [zarr[Nm1], zarr[0]], [Ed[Nm1], Es[Nm1], Eu[Nm1]], 'RK45', t_eval = np.flip(zarr))
    
    Ed = np.flip(res.y[0])
    Es = np.flip(res.y[1])
    Eu = np.flip(res.y[2])
    
    return Ed, Es, Eu 


def ocean_irradiance(hbot, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, N = 30, 
                     pt1_perc_zbot = True, use_bvp_solver = False):
    
    
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
    init_guess = .2 
    # init_guess = 0
    
    Ed1 = np.full(N, init_guess)
    Es1 = np.full(N, init_guess)
    Eu1 = np.full(N, init_guess) 
    
    ## BCs included in the initial guess arrays
    Ed1[Nm1] = Ed0 ##surface 
    Es1[Nm1] = Es0 ##surface 
    Eu1[0] = Euh ##bottom
    
    ## Default number of shots for BVP shoot method solution
    shots = 3 
    
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
        zbot_pt1perc = zbot_func(Ed0, c_wat, z_phy)
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
    
    ## Scipy solver. 
    if use_bvp_solver:
        Ed,Es,Eu = ocean_irradiance_scipy(z, Ed0, Es0, Euh, a, b, coefficients)
    
    ## Our solver.
    else :

        
        Ed=np.copy(Ed1)
        Es=np.copy(Es1)
        Eu=np.copy(Eu1)
    
        Eu0_tried = []
        Fmetric = []
        
        for jsh in range(shots) :
       # Integrate down from the top to ensure Ed(1) and Es(1) are good.
    
            if jsh == 0:
                dEu = 0 #-Eu[Nm1]
            elif jsh == 1:
            # for the first case, need some adjustment to get gradient.
                dEu = max(0.01,0.03*Es[Nm1])
                # dEu = .2
            else: 
                Jslope = (Fmetric[jsh-2]-Fmetric[jsh-1]) / (Eu0_tried[jsh-2]-Eu0_tried[jsh-1]) 
                dEu = -Fmetric[jsh-1]/Jslope
    
            dEu = max(-Eu[Nm1], dEu)
            dEu = min(1-Eu[Nm1], dEu)
            Eu[Nm1] = Eu[Nm1] + dEu
            Eu0_tried.append(Eu[Nm1])
    
            Edmid = np.zeros(Nm1)
            ## integrate Es down the water column
            ## The range does not actually go to k=-1, only to k=0. 
            ## i.e. does not include stop point of range. 
            # for k in range(Nm1-1 , -1, -1) :
            Ed, Es, Eu = Irradiance_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, 
                                        r_s, r_u, v_d, v_s, v_u, direction = 'down')
            # Ed, Es, Eu = Scipy_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, 
            #                             r_s, r_u, v_d, v_s, v_u)
                
                ## calculate a metric that indicates goodness of our shot.
                ## since Eu(bot) = 0, our metric of it is just the bottom value for Eu.
            Fmetric.append(Eu[0])
        # Eu[-1] = Eu[-1] - .0001
        # Ed, Es, Eu = Irradiance_RK4(Nm1, Ed, Es, Eu, z, a, b, c_d, b_b, b_f, 
        #                             r_s, r_u, v_d, v_s, v_u)


    return Ed, Es, Eu, z


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
    shots = 3
    
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


def ocean_irradiance_dutkiewicz(hbot, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, N = 30, 
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
    
    def E_s_z(z,zarr,c_p, c_m, E_d_z):
        """
    
        Parameters
        ----------
        z : 1-D array of length N-1
            This is the z that is chosen to represent z at in place in the column.
            It is chosen for now to be the midpoint of the layer. 
            Later could be chosen as any point but must know what layer the point exists in 
        zarr : 1-D array of length N
            The array that I have used for all the coefficients in this program 
        c_p : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
        c_m : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
    
        Returns
        -------
        E_s : 1-D array of length N-1
            The Downward diffuse irradiance solution. 
    
        """
        ##midpoint z values in general solution. 
        # E_d_z = analytical_Ed(z, c_Ed_z, Ed0) ##the downward direct irradiance at midpoints z 
        E_s = np.zeros(N-1)
        for k in range(N-1):
            E_s[k] = (c_p[k])*(np.exp((-kap_p[k])*(-(z[k] - zarr[k])))) + (c_m[k])*(r_m[k])*((np.exp((kap_m[k])*(-(z[k] - zarr[k+1]))))) + (x[k])*(E_d_z[k])
        return E_s 
    
    def E_u_z(z,zarr,c_p, c_m, E_d_z):
        """
    
        Parameters
        ----------
        z : 1-D array of length N-1
            This is the z that is chosen to represent z at in place in the column.
            It is chosen for now to be the midpoint of the layer. 
            Later could be chosen as any point but must know what layer the point exists in 
        zarr : 1-D array of length N
            The array that I have used for all the coefficients in this program 
        c_p : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
        c_m : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
    
        Returns
        -------
        E_u : 1-D array of length N-1
            The upwelling irradiance solution. 
    
        """
        ##midpoint z values in general solution. 
        # E_d_z = analytical_Ed(z, c_Ed_z, Ed0) ##the downward direct irradiance at midpoints z 
        E_u = np.zeros(N-1)
        for k in range(N-1):
            E_u[k] = (c_p[k])*(r_p[k])*(np.exp((-kap_p[k])*(-(z[k] - zarr[k])))) + (c_m[k])*(r_m[k])*((np.exp((kap_m[k])*(-(z[k] - zarr[k+1]))))) + (y[k])*(E_d_z[k])
        return E_u 
    
    
    
    ## PARAMS FROM DUTKIEWICZ 2015 
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    ##N centers
    Nm1 = N - 1  
    
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
                
    c_wat = (a_wat+b_wat)/v_d ##used for analytical 
    c_d = (a+b)/v_d
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        # zbot_pt1perc = zbot_func(Ed0, a_wat, b_wat, v_d)
        zbot_pt1perc = zbot_func(Ed0, c_wat, z_phy)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    z = Log_Trans(zbot, N) 
    z_out = np.zeros(Nm1)
    for k in range(Nm1):   
        dz = z[k+1] - z[k]  
        z_out[k] = z[k] + dz/2 
    z_out = np.flip(z_out)
    ## linear z 
    #z = np.linspace(zbot, 0, N)
    
    
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    if phy: 
        # print(z_phy)
        # print(z)
        ## FLipping the coordinates because the interpolation requires 'monotonically increasing'
        a = np.flip(np.interp(z,z_phy,a))
        b = np.flip(np.interp(z,z_phy,b))
        z = np.flip(z)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)

    b_b = .551*b 
    b_f = b - b_b 

     ## Don't know much about this const/variable 
    c = (a+b)/v_d ##used for analytical 
    c_d = c
    a2 = a[:-1]
    b2 = b[:-1]
    c_Ed_z = (a2+b2)/v_d ##used in functions below for midpoint 
    
    ##maybe it is the downward direct coefficient?
    
    ##Making the matching constant of Dutkiewicz 
    C_s = (a + r_s*b_b)/ v_s ##Cs 
    C_u = (a + r_u*b_b)/ v_u ## Cu 
    
    B_u = (r_u*b_b)/v_u 
    B_s = (r_s*b_b)/v_s 
    
    F_d = b_f / v_d  ##NOTE : Here I don't use what Dutkiewicz uses for F_d, B_d
    B_d = b_b/ v_d  ##I use the coefficient of E_d from eq. 1,2,3 of Dutkiewicz 
    
    ##Inhomogeneous solution following Dutikiewicz, Eqaution B9 
    ##first the det of M 
    det_M = 1/((c_d - C_s)*(c_d + C_u) + B_s*B_u)
    x = det_M*(-F_d*(c_d + C_u) - B_u*B_d)
    y = det_M*(-F_d*B_s + B_d*(c_d - C_s))
    
    ##now definining some stuff for the homogeneous solution 
    D = .5*(C_s + C_u + np.sqrt((C_s + C_u)**2 - 4*B_s*B_u))
    
    ##eigen values 
    kap_m = D - C_s ##kappa minus 
    kap_p = -(C_u - D )##kappa plus 
    
    ##eigen vectors 
    r_p = B_s / D ## r plus 
    r_m = B_u / D ## r minus 
    
    ##defining the exponential decays in each layer 
    ##note that because our grid has a uniform spacing 
    ##we can just let z_k+1 - z_k = dz = const. 
    ##in ROMS this might not be the case, but it is a valid choice rn 
        ##Now for making the coefficient matrix 
    ##For this problem instead of stacking the state vectors
    ##I will be interweaving them by alternating c^plus and c^minus. 
    ##This should give a tri-diagonal matrix 
    ##Which can be solved throuhg Gaussian elimination 
    
    ## A zero 2N by 2N matrix 
    A = np.zeros((2*N,2*N))
    
    for k in range(0, N-1): ##this leaves space for boudnaries at top/ bottom 
        dz = abs(z[k+1] - z[k])
        e_p = np.exp(-kap_p*dz) ##dz = const 
        e_m = np.exp(-kap_m*dz) ##dz = xonst  
    
       ##since there is only k and k+1 it is ok to start at k =0
        c1 = (e_p[k])*(1 - (r_p[k])*(r_m[k+1]))
        c2 = (r_m[k]) - (r_m[k+1])
        c3 = (1 - (r_p[k+1])*(r_m[k+1]))
        #cxy = x[k+1] - x[k]- (y[k+1] - y[k])*(r_m[k+1]) ##needed for E_d vector
        
        c4 = (1 - (r_m[k])*(r_p[k]))
        c5 = ((r_p[k+1]) - (r_p[k]))
        c6 = (e_m[k+1])*(1- (r_m[k+1])*(r_p[k])) 
        #cyx = y[k+1] - y[k] - (x[k+1] - x[k])*(r_p[k]) ##needed for E_d vector
        
        #if (k % 2) == 0: ##if odd, c1,c2,c3 will be first top of the weaved stack 
        m = 2*k + 1 ##odd numbers, k starting at zero 
        A[m,(m - 1)] = c1 
        A[m,m] = c2 
        A[m,(m+1)] = -c3
        
        n = 2*k + 2 ##even numbers k starting at zero 
        A[n,(n - 1)] = c4 
        A[n,n] = -c5 
        A[n,(n+1)] = -c6
    
    ##now some boundary conditions
    ## c_kbot^+ = 0 at bottom thus, 
    A[-1,-1] = 1
    
    
    ##at top we use Dutkiewicz et al.(2015) setting E_s0 - xE_d0 = somthing 
    A[0,0] = 1
    A[0,1] = (r_m[0])*np.exp(-(kap_m[0])*(z[1])) ## = E_s0 - x[0] * E_d0

    # E_d = analytical_Ed(z, c, Ed0)
    E_d =  np.flip(numerical_Ed(np.flip(z), np.flip(c_Ed_z), Ed0))
    
    E_d2 = np.zeros(2*N)
    
    for k in range(N-1):
        cxy = x[k+1] - x[k]- (y[k+1] - y[k])*(r_m[k+1]) ##needed for E_d vector
        cyx = y[k+1] - y[k] - (x[k+1] - x[k])*(r_p[k]) ##needed for E_d vector
        E_d2[2*k+1] = cxy*E_d[k] ##technically this is E_d at k+1
        E_d2[2*k+2] = cyx*E_d[k] ##techinically this is also
    
    ##now for setting the boundaries 
    E_d2[0] = Es0 - (x[0])*Ed0
    E_d2[-1] = 0
    
    ##solving by LU-decomp. 
    B= E_d2 
    
    lu, piv = sp.linalg.lu_factor(A)
    x_lu = sp.linalg.lu_solve((lu, piv), B)
    c_p = np.zeros(N)
    c_m = np.zeros(N)
    for i in range(2*N): 
        if (i%2) == 0: 
            c_p[int(i/2)] = x_lu[i]
        else: 
            c_m[int((i-1)/2)] = x_lu[i]
    
      
    #z_out = np.linspace(z[0] ,z[-1], N-1) ##z array for E_s_z and E_u_z 
   # z_out = np.linspace(0,z[-1] + dz/2, N-1)

    Ed = np.flip(numerical_Ed(np.flip(z_out), np.flip(c_Ed_z), Ed0))
    Es = E_s_z(z_out, z, c_p, c_m, Ed)
    Eu = E_u_z(z_out, z, c_p, c_m, Ed)

    
    # print(c_Ed_z)
    #Es = x_lu[:N]
    #Eu = x_lu[N:] 
    #Ed = analytical_Ed(zarr, c)
    
    return Ed, Es, Eu, z_out


def ocean_irradiance_dutkiewicz_ROMS(hbot, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, N = 30, 
                     pt1_perc_zbot = True):
    
    def E_s_z(Nm1, E_d_z, z, zarr, c_p, c_m):
        """
        
        Parameters
        ----------
        z : 1-D array of length N-1
            This is the z that is chosen to represent z at in place in the column.
            It is chosen for now to be the midpoint of the layer. 
            Later could be chosen as any point but must know what layer the point exists in 
        zarr : 1-D array of length N
            The array that I have used for all the coefficients in this program 
        c_p : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
        c_m : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
        
        Returns
        -------
        E_s : 1-D array of length N-1
            The Downward diffuse irradiance solution. 
        
        """
        ##midpoint z values in general solution. 
        # E_d_z = numerical_Ed(c) ##the downward direct irradiance at edges z 
        E_s = np.zeros(Nm1)
        for k in range(Nm1,0,-1):
            E_s[k-1] = (c_p[k])*(np.exp((-kap_p[k])*(-(z[k-1] - zarr[k])))) + (c_m[k])*(r_m[k])*((np.exp((kap_m[k])*(-(z[k-1] - zarr[k-1]))))) + (x[k])*(E_d_z[k])
        return E_s 


    def E_u_z(Nm1, E_d_z, z, zarr, c_p, c_m):
        """
        
        Parameters
        ----------
        z : 1-D array of length N-1
            This is the z that is chosen to represent z at in place in the column.
            It is chosen for now to be the midpoint of the layer. 
            Later could be chosen as any point but must know what layer the point exists in 
        zarr : 1-D array of length N
            The array that I have used for all the coefficients in this program 
        c_p : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
        c_m : 1-D array of length N
            The values of c^+ that I found through Gaussian Elimination of tridiagonal
        
        Returns
        -------
        E_u : 1-D array of length N-1
            The upwelling irradiance solution. 
        
        """
        ##midpoint z values in general solution. 
        # E_d_z = numerical_Ed(c) ##the downward direct irradiance at edges 
        E_u = np.zeros(Nm1)
        for k in range(Nm1,0,-1):
            
            E_u[k-1] = (c_p[k])*(r_p[k])*(np.exp((-kap_p[k])*(-(z[k-1] - zarr[k])))) + (c_m[k])*(r_m[k])*((np.exp((kap_m[k])*(-(z[k-1] - zarr[k-1]))))) + (y[k])*(E_d_z[k])
        
        return E_u 

    ## PARAMS FROM DUTKIEWICZ 2015 
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    ##N centers
    Nm1 = N - 1  
    
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
                
    c_wat = (a_wat+b_wat)/v_d ##used for analytical 
    c_d = (a+b)/v_d
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        # zbot_pt1perc = zbot_func(Ed0, a_wat, b_wat, v_d)
        zbot_pt1perc = zbot_func(Ed0, c_wat, z_phy)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    # z = Log_Trans(zbot, N) 
    ## linear z 
    zarr = np.linspace(zbot, 0, N)
    
    
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    if phy: 
        # print(z_phy)
        # print(z)
        ## FLipping the coordinates because the interpolation requires 'monotonically increasing'
        a = np.interp(zarr,z_phy,a)
        b = np.interp(zarr,z_phy,b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)

    b_b = .551*b 
    b_f = b - b_b 
    c = (a+b)/v_d ##used for analytical 
    c_d = c
    
    ##maybe it is the downward direct coefficient?
    
    ##Making the matching constant of Dutkiewicz 
    C_s = (a + r_s*b_b)/ v_s ##Cs 
    C_u = (a + r_u*b_b)/ v_u ## Cu 
    
    B_u = (r_u*b_b)/v_u 
    B_s = (r_s*b_b)/v_s 
    
    F_d = b_f / v_d  ##NOTE : Here I don't use what Dutkiewicz uses for F_d, B_d
    B_d = b_b/ v_d  ##I use the coefficient of E_d from eq. 1,2,3 of Dutkiewicz 
    
    ##Inhomogeneous solution following Dutikiewicz, Eqaution B9 
    ##first the det of M 
    det_M = 1/((c_d - C_s)*(c_d + C_u) + B_s*B_u)
    x = det_M*(-F_d*(c_d + C_u) - B_u*B_d)
    y = det_M*(-F_d*B_s + B_d*(c_d - C_s))
    
    ##now definining some stuff for the homogeneous solution 
    D = .5*(C_s + C_u + np.sqrt((C_s + C_u)**2 - 4*B_s*B_u))
    
    ##eigen values 
    kap_m = D - C_s ##kappa minus 
    kap_p = -(C_u - D )##kappa plus 
    
    ##eigen vectors 
    r_p = B_s / D ## r plus 
    r_m = B_u / D ## r minus 
    
    ##defining the exponential decays in each layer 
    ##note that because our grid has a uniform spacing 
    ##we can just let z_k+1 - z_k = dz = const. 
    ##in ROMS this might not be the case, but it is a valid choice rn 
    dz = abs(zarr[1] - zarr[0])
    e_p = np.exp(-kap_p*dz) ##dz = const 
    e_m = np.exp(-kap_m*dz) ##dz = xonst  
    
    ## A zero 2N by 2N matrix 
    A = np.zeros((2*N,2*N))
    
    for k in range(Nm1,0,-1): ##this leaves space for boudnaries at top/ bottom 
        ##since there is only k and k+1 it is ok to start at k =0
        c1 = (e_p[k])*(1 - (r_p[k])*(r_m[k-1]))
        c2 = (r_m[k]) - (r_m[k-1])
        c3 = (1 - (r_p[k-1])*(r_m[k-1]))
        #cxy = x[k+1] - x[k]- (y[k+1] - y[k])*(r_m[k+1]) ##needed for E_d vector
        
        c4 = (1 - (r_m[k])*(r_p[k]))
        c5 = ((r_p[k-1]) - (r_p[k]))
        c6 = (e_m[k-1])*(1- (r_m[k-1])*(r_p[k])) 
        #cyx = y[k+1] - y[k] - (x[k+1] - x[k])*(r_p[k]) ##needed for E_d vector
        
        #if (k % 2) == 0: ##if odd, c1,c2,c3 will be first top of the weaved stack 
        m = 2*k  ##odd numbers, k starting at zero 
        A[m,(m + 1)] = c1 
        A[m,m] = c2 
        A[m,(m - 1)] = -c3
        
        n = 2*k - 1 ##even numbers k starting at zero 
        A[n,(n + 1)] = c4 
        A[n,n] = -c5 
        A[n,(n - 1)] = -c6
    
    ##now some boundary conditions
    ## c_kbot^+ = 0 at bottom thus, 
    A[0,0] = 1
    
    
    ##at top we use Dutkiewicz et al.(2015) setting E_s0 - xE_d0 = somthing 
    A[-1,-1] = 1
    A[-1,-2] = (r_m[-1])*np.exp(-(kap_m[-1])*(zarr[-2])) ## = E_s0 - x[0] * E_d0
    
    
    ##YAY A is done 
    
    ##no for E_d 
    
    E_d = numerical_Ed(zarr, c_d, Ed0)
    
    E_d2 = np.zeros(2*N)
    
    for k in range(Nm1,0,-1):
        cxy = x[k-1] - x[k]- (y[k-1] - y[k])*(r_m[k-1]) ##needed for E_d vector
        cyx = y[k-1] - y[k] - (x[k-1] - x[k])*(r_p[k]) ##needed for E_d vector
        E_d2[2*k] = cxy*E_d[k] ##technically this is E_d at k+1
        E_d2[2*k-1] = cyx*E_d[k] ##techinically this is also
    
    ##now for setting the boundaries 
    E_d2[-1] = Es0 - (x[-1])*Ed0
    E_d2[0] = 0
    
    
    
    ##solving by LU-decomp. 
    B = E_d2 
    
    lu, piv = sp.linalg.lu_factor(A)
    x_lu = sp.linalg.lu_solve((lu, piv), B)
    
    c_p = np.zeros(N)
    c_m = np.zeros(N)
    # for i in range(2*Nm1,0,-1): 
    #     if (i%2) == 0: 
    #         c_m[int(i/2)] = x_lu[i]
    #     else: 
    #         c_p[int((i+1)/2)] = x_lu[i]
    for i in range(N,0,-1):
          c_p[i-1] = x_lu[2*(i)-1]
          c_m[i-1] = x_lu[2*(i-1)]
    

    
 
    
    z = np.linspace(zarr[Nm1], 0, Nm1) ##z array for E_s_z and E_u_z 
    
    Ed_zarr = numerical_Ed(zarr, c_d, Ed0)
    Es = E_s_z(Nm1, Ed_zarr, z, zarr, c_p, c_m)
    Eu = E_u_z(Nm1, Ed_zarr, z, zarr, c_p, c_m)
    Ed = numerical_Ed(z, c_d, Ed0)
    
    
    return Ed, Es, Eu, z
    

def artificial_phy_prof(z,loc,width,conc):

    # prof = conc*(1 + np.tanh((z-loc)/width)) 
    
    prof = conc* np.exp(-((z-loc)**2)/(2*width**2))
    
    return prof


def Demo(): 
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

    phy_prof = artificial_phy_prof(z, -90, 40,1.5)
    # ROMS_point = np.genfromtxt('ChrisData_good_point.csv', delimiter=',')
    # phy_prof = ROMS_point[1:,2]
    # print(phy_prof)
    # z = ROMS_point[1:,0]
    # phy_prof = np.full(len(z), 1)
    
    ab_wat = abscat(lam, 'water')
    
    a_phy, b_phy = abscat(lam, 'Syn')
    
    # a_phy = .01
    # b_phy =  .0039 
    
    ## Define the Phytoplankton class.
    phy = Phy(z, phy_prof, a_phy, b_phy)

    # Ed0 = .7
    # Es0 = 1 - Ed0
    # Euh = 0 

    zbot = z[0]

    # Ed, Es, Eu, zarr = ocean_irradiance_dutkiewicz(zbot,PI.Ed0,PI.Es0,PI.Euh,ab_wat, PI.coefficients, 
    #                                                 phy=phy, N=N, pt1_perc_zbot = False)
    
    # dev = .001
    # zbot = zarr[Es > dev][-1]
    
    # Ed, Es, Eu, zarr, c_Ed_z = ocean_irradiance_dutkiewicz(zbot,PI.Ed0,PI.Es0,PI.Euh,ab_wat, PI.coefficients, 
    #                                             phy=phy, N=N, pt1_perc_zbot = False)
                                        
    Ed, Es, Eu, zarr = ocean_irradiance(zbot,PI.Ed0,PI.Es0,PI.Euh,ab_wat, PI.coefficients, 
                                        phy=phy, N=N, pt1_perc_zbot = False,
                                        use_bvp_solver = False)

    ## Plotting the Results
    #-------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    
    ax1.plot(phy_prof, z)
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
    return zarr, Ed


#-------------------------------------MAIN-------------------------------------

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Ocean Irradiance Fucntion')

    parser.add_argument('--demo', action='store_true',
                        help='Run the demo and plot')
    
    args = parser.parse_args()
    
    # if args.demo: 
    zarr, Ed = Demo()
    
    #Ed_redo = np.flip(numerical_Ed(np.flip(zarr), np.flip(c_Ed_z), .7))





###################################ARCHIVED WORK##############################

# def analytical_Ed(zarr,c,E_d_0):
    
#     Np1 = len(zarr)  # N+1 edges
#     N=Np1-1          # N cell centers on which c is defined
#     Nm1 = N - 1
   
#     if len(c)!=N:
#         print('incompatible lengths of c and zarr')

#     Ed = np.zeros(Np1)
    
#     Ed[N]=E_d_0
#     for k in range(Nm1,-1,-1): ##k goes from Nm1 --> 0
#         Ed[k] = Ed[k+1]*np.exp(-c[k]*(zarr[k+1]-zarr[k]))

#     return Ed
