"""
Created on Thu Dec 10 08:52:28 2020

@author: Miles Miller
"""
import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt
from scipy import integrate



class Phy:
    """
    The Phytoplankton Class includes:
        --> z coordinates corresponding to phytoplankton concentrations [m (< 0)]
                -- z should be a column vector
        --> Phytoplankton chlorophyll concentrations [mg chl-a m^-3]
                --2-D array with columns corresponding to diff. phytoplankton
        --> Equivalent Spherical Diameter
                -- 1-D Array corresponding to the number of phytoplankton
        --> absorbtion coefficient, a [m^-1]
                -- Vector corresponding to different phytoplankton
        --> backscatter coefficient, b [m^-1]
                -- Vector corresponding to different phytoplankton
    """
    def __init__(self, z, phy, esd, a, b):
        self.z = z
        self.phy = phy
        self.esd = esd
        self.a = a
        self.b = b
        self.Nz = len(z)
        if phy.ndim == 1:
            self.Nphy = 1
        elif phy.ndim == 2:
            self.Nphy = np.shape(phy)[1]
            assert self.Nphy == len(a) and self.Nphy == len(b)
            assert self.Nphy == len(esd)

        assert np.shape(phy)[0] == self.Nz
        
        return 


class CDOM_sal:
    """
    The CDOM_sal Class includes:
        --> z coordinates corresponding to concentrations [m (< 0)]
                -- z should be a column vector
        --> salinity concentrations []
                --1-D array same size as z
        --> given wavelength
                -- float 
    The CDOM is calculated from the salinity following equations 2 and 3 from 
    "Remote-sensing reflectance and true colour produced by a coupled
     hydrodynamic, optical, sediment, biogeochemical model of the 
     Great Barrier Reef, Australia: Comparison with satellite data"
    by Mark E. Baird 2016

    """
    def __init__(self, z, salt, wavelength):
        self.z = z
        self.salt = salt
        self.wavelength = wavelength

        ## calculate cdom from salitnity following Baird 2016 
        ## Edited by Jonathan such that absorbtion = 0 at 34.
        salt[salt >= 34] = 34
        a_443 = -0.0332*salt + 1.1288
        self.a = a_443*np.exp(-0.012*(wavelength - 443))
        
        return 


class CDOM_dens: 
    """
    The CDOM_dens class includes: 
        --> z coordinates coorespnding to ceoncentration levels [m (<0)]. 
        --> CDOM density
        --> Absorption coefficient in units of m^2(mmolC)^-1.  
        --> CDOM to Carbon ratio. This would be CDOM [C] = CDOM [mg/m^3] * ratio. 
    This cdom uses observed CDOM concentration to calculate the absorption due to 
    CDOM. 

    """

    def __init__(self, z, cdom, CDOM2C, wavelength): 
        self.z = z
        self.cdom = cdom 
        self.CDOM2C = CDOM2C
        self.wavelength = wavelength

        ## [This absorption curve is given by Dutkiweicz et al. 2015 equation 15. ] 
        ## [Value given is in table 1.]
        c_cdom = 0.18 ## [m^2 (mmol C)^-1]
        lam_0 = 450 ## [nm]
        s_cdom = 0.021 ## [(nm)^-1]
        ## [Equation 15.]
        self.a = c_cdom * np.exp(-s_cdom*(wavelength - lam_0)) 
        
        

        return 
        

class CDOM_refa: 
    """
    This class is meant for the approximation of constant CDOM absorption. This takes a reference absorption and calculates the 
    absorption spectrum for CDOM.

    Parameters
    ----------
    z: 1-D Array [N]
        The z coordinates on which the absorption due to cdom is defined.
    cdom_refa: 1-D Array [N]
        The reference absorption value for the CDOM absorption spoectrum. Choose something small like 400-450nm.
    lam0: Float
        The reference absorption corresponding wavelength value. 
    wavelength: Float
        The wavelength of the desired absorption. 
    fraca: Optional, Float
        Default is 1.0. This is the fraction of the reference absorption to be used. For example if 
        50% of the absorption at the reference wavelength is due to CDOM then fraca would equal 0.5.
    
    Returns
    -------
    a: 1-D Array [N]
        An array filled with the scalar constant CDOM absroption at the provided wavelength. 

    """

    def __init__(self, z, cdom_refa, lam0, wavelength, fraca=1.0): 

        self.z = z
        self.cdom_refa = cdom_refa
        self.lam0 = lam0 
        self.wavelength = wavelength
        self.fraca = fraca

        ## [This absorption curve is given by Dutkiweicz et al. 2015 equation 15. ] 
        ## [Value given is in table 1.]
        s_cdom = 0.021 ## [(nm)^-1]
        ## [Equation 15.]
        a = fraca *cdom_refa * np.exp(-s_cdom*(wavelength - lam0)) 
        self.a = a
        
        return 


class Det: 
    """
    This class is the class for the detritus. It includes its absorption, 
    scattering, and backscatter coefficients. It also includes the
    concentration profile of detritus and the associated depths.
    """

    def __init__(self, z, det, a, b, b_b): 
        self.z = z
        self.det = det
        self.a = a
        self.b = b 
        self.b_b = b_b

        return 


def Backscatter_Ratio(esd):
    """ 
    This calculates the backscatter ratio for a given equivalent spherical diameter. 
    This equation comes from equation 6 of the paper:

    "Spectral backscattering properties of marine phytoplankton cultures"
    by 
    Amanda L. Whitmire, W. Scott Pegau, Lee Karp-Boss, Emmanuel Boss, and Timothy J. Cowles1
    2010

    Parameters
    ----------
    esd: Float
        The equivalent spherical diameter of the given phytoplankton species. Should be 
        in micro meters.
    
    Returns
    -------
    bb_r: Float
        The backscatter ratio.
    """
   
    bb_r = (4.390*10**-3)*esd**0.432
    
    return bb_r


def Calc_Abscat_Grid(hbot, ab_wat, N, phy=None, CDOM_refa=None, det=None, grid='log', pt1_perc_zbot=True, pt1_perc_phy=True):

    """
    This function calculates the irradiance grid, calculates the absorption, scattering, and backscatter 
    profiles and then interpolates those onto the irradiance grid. 

    It assumes that the vertical grids for the different consituents are the same. 

    Parameters
    ----------

    Returns
    -------
    """

    ##unpacking the ab_wat_tuple 
    a_wat,b_wat = ab_wat 
    a = a_wat
    b = b_wat 
    
    ## backscattering due to phy
    b_b_phy = 0
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
        
        ## equivalent spherical diameter
        esd = phy.esd
        
        ## Just one phytoplankton species
        if Nphy == 1 : 
            ## The back scatter ratio
            bb_r = OI.Backscatter_Ratio(esd)    
            a = a + phy_prof * a_phy
            b = b + phy_prof * b_phy
            b_b_phy = b_b_phy + phy_prof * b_phy * bb_r
            
        ## More than one species
        elif Nphy > 1 : 
            for k in range(Nphy):
                ## The back scatter ratio
                bb_r = Backscatter_Ratio(esd[k])    
                a = a + phy_prof[:,k] * a_phy[k]  
                b = b + phy_prof[:,k] * b_phy[k]
                b_b_phy = b_b_phy + phy_prof[:,k] * b_phy[k] * bb_r

    ## [Inclusion of CDOM via its reference absorption.]
    if CDOM_refa: 
        ## [Unpack.]
        z_cdom = CDOM_refa.z
        a_cdom = CDOM_refa.a

        ## [Assert that the z grid is the same for cdom and phy.]
        assert ~np.any(z_phy != z_cdom), "Please ensure that consituents are on the same vertical grid."

        ## [Including CDOM in the absorption.]
        a += a_cdom

    ## [The inclusion of detritus.]
    if det: 
        ## [unpack the onject.]
        z_det = det.z
        det_prof = det.det
        a_det = det.a
        b_det = det.b
        b_b_det = det.b_b


        ## [Assert that the z grid is the same for cdom and phy.]
        assert ~np.any(z_phy != z_det), "Please ensure that consituents are on the same vertical grid."
        
        ## [Inclusion of detritus.]
        a += a_det * det_prof
        b += b_det * det_prof
        b_b_det = b_b_det * det_prof
            

    ## [Forming the irradiance grid.]
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        if pt1_perc_phy == True:
            zbot_pt1perc = OI.zbot_func(Ed0, a, b, v_d, phy=True, z=z_phy) 
        else:    
            zbot_pt1perc = OI.zbot_func(Ed0, a_wat, b_wat, v_d, phy=False)
        if zbot_pt1perc == None:
            print('bad pt1 perc light level')
            zbot_pt1perc = -100
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 

    ## [construct the vertical grid to the designated depth.]
    if grid == 'log': 
        ## log transformed z grid.
        z = OI.Log_Trans(zbot, N) 
    elif grid == 'linear'
        ## linear z 
        z = np.linspace(zbot, 0, N)
    
    ## Interpolating a,b vectors from z_phy to z.
    if phy: 
        a = np.interp(z,z_phy,a)
        b = np.interp(z,z_phy,b)
        b_b_phy = np.interp(z, z_phy, b_b_phy)
        if det: 
            b_b_det = np.interp(z, z_phy, b_b_det)
    elif CDOM_sal or CDOM_dens or CDOM_refa: 
        a = np.interp(z,z_cdom,a)
        ## no scattering for cdom, thus just water scattering.
        b = np.full(N, b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        
    b_b_wat = .551*b_wat

    b_b = b_b_wat + b_b_phy

    return z, a, b, b_b
 

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


def zbot_func(Ed0, a, b, v_d, light_frac = .01, phy=False, z=None):
    """
    Finds the zbot for at which light ha attenuated to .1% of its surface value 
    for water only coeffients

    Parameters
    ----------
    Ed0 : Float 
        Initial Value for E_d. 
    a: 1-D Array
        The total absorbtion necessary for the calculation of Ed.
    b: 1-D Array
        The total scattering necessary for the calculation of Ed.
    v_d: Float
        The averaging of the cosine in the Ed equation. 
    light_frac: Float
        Optional, the default is .01. 
    phy: Boolean
        Optional, the default is False.
    z: 1-D Array
        The z-coordinates corresponding to the absorbtion and scattering arrays
        in the case that they are not constants. 
    Returns
    -------
    zbot : Float 
        .01% light level zbot. 

    """
    
    ## For now let the scattering be zero. 
    #b = 0
    c = (a+b) / v_d
    #zbots = np.linspace(-1000, 0, 10000) 
    zbots = Log_Trans(-1500, 5000) 
    
    if phy==True: 
        c = np.interp(zbots, z, c)
        Ed = numerical_Ed(zbots, c, Ed0)
    else:
        Ed = analytical_Ed(zbots, c, Ed0)
    ## The flipping is so the iteration starts at the surface.
    for k, Ed_i in enumerate(np.flip(Ed)) :
        EdoEd0 = Ed_i / Ed0
        #print(EdoEd0, light_frac, np.flip(zbots)[k], phy, Ed_i)
        if EdoEd0 < light_frac :
            zbot = np.flip(zbots)[k] 
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
    zarr_exp = np.linspace(np.log(abs(zbot-1)),0, Nlayers)
    zarr = 1-np.exp(zarr_exp)
    
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
    z_mesh = res.x

    Ed = np.interp(zarr, z_mesh, y[0])
    Es = np.interp(zarr, z_mesh, y[1])
    Eu = np.interp(zarr, z_mesh, y[2])
    
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
            
            ## Setting very small Eu to zero. 
            if Eu[k+1] < 1e-10:
                Eu[k+1] = 0
            

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
        zbot_pt1perc = zbot_func(Ed0, c_wat)
        #print(zbot_pt1perc)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    z = Log_Trans(zbot, N) 
    ## linear z 
    #z = np.linspace(zbot, 0, N)
    
    
    
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


def ocean_irradiance_shoot_up(hbot, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, CDOM=None, N = 30, pt1_perc_zbot = True, pt1_perc_phy = True):
    
    
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
    pt1_perc_phy: Boolean, default is True
        The .1% light level is used to set zbot, but this flag calculates the 
        .1% light level with phytoplankton. The pt1_perc_zbot flag must also be True for this
        to work.

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
    
    ## backscattering due to phy
    b_b_phy = 0
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
        
        ## equivalent spherical diameter
        esd = phy.esd
        ## Just one phytoplankton species
        if Nphy == 1 : 
            ## The back scatter ratio
            bb_r = Backscatter_Ratio(esd)    
            a = a + phy_prof * a_phy
            b = b + phy_prof * b_phy
            b_b_phy = b_b_phy + phy_prof * b_phy * bb_r
            
        ## More than one species
        elif Nphy > 1 : 
            for k in range(Nphy):
                ## The back scatter ratio
                bb_r = Backscatter_Ratio(esd[k])    
                a = a + phy_prof[:,k] * a_phy[k]  
                b = b + phy_prof[:,k] * b_phy[k]
                b_b_phy = b_b_phy + phy_prof[:,k] * b_phy[k] * bb_r

    ## Inclusion of CDOM
    if CDOM: 
        ## unpacking the object
        ## For now it is assumed that phy and cdom share same grid.
        if phy:
            a_cdom = CDOM.a
            a = a + a_cdom 
        ## This is for only CDOM no phy
        else:
            a_cdom = CDOM.a
            z_cdom = CDOM.z 
            a = a + a_cdom 

    
    ## Irradiance Grid Stuff
    ## If pt1_perc_zbot is True
    if pt1_perc_zbot == True :
        ## Finding the zbot at the .1% light level. 
        if pt1_perc_phy == True:
            c_d = (a+b)/v_d
            zbot_pt1perc = zbot_func(Ed0, a, b, v_d, phy=True, z=z_phy) 
        else:    
            c_wat = (a_wat + b_wat)/v_d
            zbot_pt1perc = zbot_func(Ed0, a_wat, b_wat, v_d, phy=False)
        if zbot_pt1perc == None:
            print('bad pt1 perc light level')
            zbot_pt1perc = -100
        #print(zbot_pt1perc)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    z = Log_Trans(zbot, N) 
    ## linear z 
    #z = np.linspace(zbot, 0, N)
    
    
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    if phy: 
        a = np.interp(z,z_phy,a)
        b = np.interp(z,z_phy,b)
        b_b_phy = np.interp(z, z_phy, b_b_phy)
    elif CDOM: 
        a = np.interp(z,z_cdom,a)
        ## no scattering for cdom, thus just water scattering.
        b = np.full(N, b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        
    ##coefficient of downward direct irradiance 
    c_d = (a+b)/v_d 
    # if N != len(Ed1) or N !=len(a)+1 or N !=len(b)+1 :
    #     print('lengths of z and Ed must be the same, and a&b should be 1 less.')
        

    b_b_wat = .551*b_wat

    b_b = b_b_wat + b_b_phy
    b_f = b - b_b 
    
    #print('a:', a, 'b:', b)
    #print('b_f:', b_f, 'b_b:', b_b)

        
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


def ocean_irradiance_shoot_fp(hbot, Ed0, Es0, Euh, ab_wat, coefficients, phy = None, N = 30, 
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
    
    
    def get_fp(Ed_val,z,c_d,Ed0):
        """
        This function gets the fixed point index as a function of the desired 
        value of Ed at which the fixed poiunt should be taken.
    
        """
        
        Ed = numerical_Ed(z, c_d, Ed0)
        for k, Edi in enumerate(Ed):
            if Edi >= Ed_val: 
                fpi = k 
                break
        
        return fpi 
    
    
    
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
        zbot_pt1perc = zbot_func(Ed0, c_wat)
        print(zbot_pt1perc)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    z = Log_Trans(zbot, N) 
    ## linear z 
    # z = np.linspace(zbot, 0, N)

    
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
    
    Ed_val = .1
    fpi = get_fp(Ed_val,z,c_d,Ed0)
    fp = z[fpi]
        

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
        zbot_pt1perc = zbot_func(Ed0, c_wat)
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
        zbot_pt1perc = zbot_func(Ed0, c_wat)
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
    

def artificial_phy_prof(z,loc,width,conc, prof_type='gauss'):

    if prof_type == 'tan': 
        prof= conc*(1 + np.tanh((z-loc)/width)) 
    elif prof_type == 'gauss': 
        prof = conc* np.exp(-((z-loc)**2)/(2*width**2))
     
    return prof


def Demo(method='shoot_up'): 
    ## Demonstration of finding the irradiance profiles for an artificial 
    ## phytoplankton concentration profile
    
    ## Absorbtion and scattering coefficients are taken from Dutkiwicz 2015 
    ## for the Diatom species. 
    import matplotlib.pyplot as plt
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    PI = Param_Init()
    
    N = 200
    Nm1 = N-1 
    lams = [443, 410]
    
    z = np.linspace(-1000,0,N)

    phy_prof = artificial_phy_prof(z, -10, 20, 1, prof_type = 'gauss')
    #phy_prof = np.full(200, 1) 
    #phy_prof = np.append(phy_prof1, np.full(100, 1))
    # ROMS_point = np.genfromtxt('ChrisData_good_point.csv', delimiter=',')
    # phy_prof = ROMS_point[1:,2]
    # print(phy_prof)
    # z = ROMS_point[1:,0]
    # phy_prof = np.full(len(z), 1)
    fig, axes = plt.subplots(1, 3, sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    
    for k, lam in enumerate(lams):
        ab_wat = abscat(lam, 'water')
    
        a_phy, b_phy = abscat(lam, 'Syn')
        
        
        ## Define the Phytoplankton class.
        phy = Phy(z, phy_prof, .02,a_phy, b_phy)
    
        ## Salinity for CDOM
        salt = np.full(N, 34.5)
        ## define the CDOM class
        cdom = CDOM(z,salt,lam) 
        
    
        ## The fixed point position: 
        fp = -50
        fpi =0
    
        zbot = z[0]
        
        if method == 'shoot_up':
            
            Ed, Es, Eu, zarr = ocean_irradiance_shoot_up(zbot,PI.Ed0,PI.Es0,PI.Euh,ab_wat, PI.coefficients, 
                                                phy=phy, CDOM=None, N=N, pt1_perc_zbot = True, pt1_perc_phy=False)
         
        if method == 'shoot_down':
            Ed, Es, Eu, zarr = ocean_irradiance(zbot, PI.Ed0, PI.Es0, PI.Euh, ab_wat,  PI.coefficients,
                                                   phy=phy, N=N, pt1_perc_zbot=True)
            
        if method == 'shoot_fp': 
            Ed, Es, Eu, zarr, fpi = ocean_irradiance_shoot_fp(zbot, fp, fpi, PI.Ed0, PI.Es0, PI.Euh, 
                                                         ab_wat, PI.coefficients, phy = phy, N = N, 
                                                         pt1_perc_zbot = True)
            
        if method == 'scipy':
            Ed, Es, Eu, zarr = ocean_irradiance(zbot, PI.Ed0, PI.Es0, PI.Euh, ab_wat,  PI.coefficients,
                                                   phy=phy, N=N, pt1_perc_zbot=False, use_bvp_solver=True)    
    
        if method == 'dut':
            Ed, Es, Eu, zarr = ocean_irradiance_dutkiewicz(zbot, PI.Ed0, PI.Es0, PI.Euh, ab_wat,  
                                                              PI.coefficients, phy=phy, N=N, pt1_perc_zbot=False)
        ## Plotting the Results
        #-------------------------------------------------------------------------
        
        markers = ['-', '-'] 
        Ed_c = 'g'
        Es_c = 'b'
        Eu_c = 'r'
        if lam == 410:
            ax2.plot(Ed, zarr, label=f'Ed', color = Ed_c, ls=markers[k] )
            ax2.plot(Es, zarr, label=f'Es', color = Es_c, ls=markers[k])
            ax2.plot(Eu, zarr, label=f'Eu', color = Eu_c, ls= markers[k])
        if lam == 551:
            ax3.plot(Ed, zarr, label=f'Ed', color = Ed_c, ls=markers[k] )
            ax3.plot(Es, zarr, label=f'Es', color = Es_c, ls=markers[k])
            ax3.plot(Eu, zarr, label=f'Eu', color = Eu_c, ls= markers[k])

    ax1.plot(phy_prof, z)
    if method == 'shoot_fp':
        ax1.hlines( zarr[fpi], min(phy_prof),max(phy_prof), color='r')
    ax1.set_xlabel('Concentration [mg Chl-a m^-3]')
    ax1.set_ylabel('Z [m]')
    ax1.set_title('Phytoplankton Concentration')
    ax1.grid()
 
    ax2.set_xlim(-0.1,1)
    ax2.set_xlabel('Irradiance')
    ax2.set_title(r'Irradiance Profiles $\lambda=443$')
    ax2.legend()
    ax2.grid()

    ax3.set_xlim(-0.1,1)
    ax3.set_xlabel('Irradiance')
    ax3.set_title(r'Irradiance Profiles $\lambda=551$')
    ax3.legend()
    ax3.grid()
    
    fig.show()
    
    
    return zarr, Ed, Es, Eu


#-------------------------------------MAIN-------------------------------------

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Ocean Irradiance Fucntion')

    parser.add_argument('--demo', action='store_true',
                        help='Run the demo and plot')
    
    args = parser.parse_args()
    
    if args.demo: 
        zarr, Ed, Es, Eu = Demo('shoot_up')
    
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