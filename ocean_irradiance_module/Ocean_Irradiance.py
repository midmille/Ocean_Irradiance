"""
Created on Thu Dec 10 08:52:28 2020

@author: Miles Miller
"""
import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt
from scipy import integrate
import ocean_irradiance_module.Wavelength_To_RGB as W2RGB



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


class CDOM_chla: 
    """
    This class is meant for the approximation of constant CDOM absorption. This takes a reference absorption and calculates the 
    absorption spectrum for CDOM.

    This calculates the cdom concentration as a function of chla concentration following: 

    Effect of a nonuniform vertical profile of chlorophylla concentration on remote-sensing reflectance of the ocean

    author: Malgorzata Stramska an dDariusz Stramski

    Applied Optics Vol. 44, Issue 9 pp 1735-1747 2005 
    """

    def __init__(self, z, chla, wavelength): 

        self.z = z
        self.chla = chla
        self.wavelength = wavelength
        
        self.a = 0.012*(chla**0.65) *np.exp(-0.014*(wavelength - 440))

        
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

def Backscatter_Ratio_2(chla): 
    """
    This function comes from 
    Bidirectional reflectances of ocean waters: accounting for Raman emission and varying particle
    scattering phase function
    
    Authors: Andre Morel, David Antoine, and Bernard Gentili

    Applied Optics Val 41 Issue 30 
    """
    
    bb_r = 0.002 + (0.01*(0.5 - 0.25* np.log10(chla)))

    return bb_r


def Calc_Abscat_Grid(hbot, ab_wat, N, Ed0, coefficients, phy=None, CDOM_refa=None, CDOM_chla=None, det=None, grid='log', pt1_perc_zbot=True, pt1_perc_phy=True):

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
    b_b_wat = 0.551*b_wat


    r_s, r_u, v_d, v_s, v_u = coefficients
    
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
#            bb_r = Backscatter_Ratio(esd)
            bb_r = Backscatter_Ratio_2(phy_prof)
            a = a + phy_prof * a_phy
            b = b + phy_prof * b_phy
            b_b_phy = b_b_phy + phy_prof * b_phy * bb_r
            
        ## More than one species
        elif Nphy > 1 : 
            
            bb_r = Backscatter_Ratio_2(np.sum(phy_prof[:,:], axis=1))
            for k in range(Nphy):
                ## The back scatter ratio
#                bb_r = Backscatter_Ratio(esd[k])    
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

    ## [Inclusion of CDOM via chla concen.]
    if CDOM_chla: 
        ## [Unpack.]
        z_cdom = CDOM_chla.z
        a_cdom = CDOM_chla.a

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
        
        light_frac_phy = 0.001
        light_frac_wat = 0.001
        ## Finding the zbot at the .1% light level. 
        if pt1_perc_phy == True:
            zbot_pt1perc = zbot_func(hbot, Ed0, a, b_b_phy, coefficients, light_frac=light_frac_phy, phy=True, z=z_phy)
        else:    
            zbot_pt1perc = zbot_func(hbot, Ed0, a_wat, b_b_wat, coefficients, light_frac=light_frac_wat, phy=False)
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
        z = Log_Trans_Grid(zbot, N) 
    elif grid == 'linear': 
        ## linear z 
        z = np.linspace(zbot, 0, N)
    else: 
        raise ValueError("Invalid grid keyword, either 'log' or 'linear'")
    
    ## Interpolating a,b vectors from z_phy to z.
    if phy: 
        a = np.interp(z,z_phy,a)
        b = np.interp(z,z_phy,b)
        b_b_phy = np.interp(z, z_phy, b_b_phy)
        if det: 
            b_b_det = np.interp(z, z_phy, b_b_det)
    elif CDOM_refa: 
        a = np.interp(z,z_cdom,a)
        ## no scattering for cdom, thus just water scattering.
        b = np.full(N, b)
    else: 
        a = np.full(N, a)
        b = np.full(N, b)
        

    ## [The total back scattering.]
    b_b = b_b_wat + b_b_phy

    ## [The total forward scattering component.]
    b_f = b - b_b 

    return z, a, b, b_b, b_f
 

def analytical_Ed_3stream(zarr, Ed0, a, b, coefficients): 
    
    
    """
    Parameters
    ----------
    zarr : vertical 1-D array
        The vertical array of the depth from 0 --> negative depth.

    Returns
    -------
    Downward Direct Irradiance
    

    """
    r_s, r_u, v_d, v_s, v_u = coefficients
   
    Ed = Ed0*np.exp(((a+b)/v_d)*zarr)
    return Ed 


def analytical_Ed_2stream(zarr, Ed0, a, b_b, coefficients): 
    
    
    """
    Parameters
    ----------
    zarr : vertical 1-D array
        The vertical array of the depth from 0 --> negative depth.

    Returns
    -------
    Downward Direct Irradiance
    

    """
    r_s, r_u, v_d, v_s, v_u = coefficients
   
    Ed = Ed0*np.exp(((a+b_b)/v_d)*zarr)
    return Ed 


def numerical_Ed_3stream(z, Ed0, a, b, coefficients):
    
    N = len(z)
    Nm1 = N - 1
    
    Ed = np.zeros(N)
    Ed[Nm1] = Ed0

    ## [The function for the rhs of 

    for k in range(Nm1, 0, -1):
        
        dz = z[k-1] - z[k]
        
        Ed[k-1] = RK45(dEddz_3stream, dz, z[k], Ed[k], a[k], b[k], None, None, coefficients)

    return Ed


def numerical_Ed_2stream(z, Ed0, a, b_b, coefficients):
    
    N = len(z)
    Nm1 = N - 1
    
    Ed = np.zeros(N)
    Ed[Nm1] = Ed0

    ## [The function for the rhs of 

    for k in range(Nm1, 0, -1):
        
        dz = z[k-1] - z[k]
        
        Ed[k-1] = RK45(dEddz_2stream, dz, z[k], Ed[k], a[k], None, b_b[k], None, coefficients)

    return Ed


def zbot_func(hbot, Ed0, a, b_b, coefficients, light_frac = .01, phy=False, z=None):
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
    zlim = -1000

    zbots = Log_Trans_Grid(max(zlim,hbot), 5000) 
    
    if phy==True: 
        a = np.interp(zbots, z, a)
        b_b = np.interp(zbots, z, b_b)
        Ed = numerical_Ed_2stream(zbots, Ed0, a, b_b, coefficients)
    else:
        Ed = analytical_Ed_2stream(zbots, Ed0, a, b_b, coefficients)
    ## The flipping is so the iteration starts at the surface.
    for k, Ed_i in enumerate(np.flip(Ed)) :
        EdoEd0 = Ed_i / Ed0
        if EdoEd0 < light_frac :
            zbot = np.flip(zbots)[k] 
            return zbot
    if zlim<hbot: 
        return hbot-10
        
    return 
   
        
def Log_Trans_Grid(zbot,Nlayers):
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


def RK45(dEdz, dz, z, E, a, b, b_b, b_f, coefficients): 
    """

    This is the RK45 method clasical
    """

    k1 = dEdz(z, E, a, b, b_b, b_f, coefficients) 
    k2 = dEdz(z+dz/2, E + dz*(k1/2), a, b, b_b, b_f, coefficients) 
    k3 = dEdz(z+dz/2, E + dz*(k2/2), a, b, b_b, b_f, coefficients) 
    k4 = dEdz(z+dz, E + dz*k3, a, b, b_b, b_f, coefficients) 
    
    Ep1 = E + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dz 

    return Ep1


def dEddz_3stream(z, Ed, a, b, b_b, b_f, coefficients): 
    """

    This the rhs of the ODE for the independent downward direct stream for the three stream
    method. 

    """

    r_s, r_u, v_d, v_s, v_u = coefficients

    dEddz = ((a+b)/v_d) * Ed

    return dEddz


def dEddz_2stream(z, Ed, a, b, b_b, b_f, coefficients): 
    """

    This the rhs of the ODE for the independent downward direct stream for the three stream
    method. 

    """

    r_s, r_u, v_d, v_s, v_u = coefficients

    dEddz = ((a+b_b)/v_d) * Ed

    return dEddz



def dEdz_3stream(z, E, a, b, b_b, b_f, coefficients): 
    """
    The rhs of the derivative of the equations for the irradiance model. This three stream model comes
    from Dutkiewicz et. al (2015) and, unlike Dutkiewicz, assumes a negative z grid.
    
    Parameters
    ----------
    E: 1-D Array, [2]
        E[0] = Es, E[1] = Eu

    Returns 
    -------
    dEdz: 1-D Array, [2]
        The solution for the rhs derivatives of the irradiance system.
    """

    r_s, r_u, v_d, v_s, v_u = coefficients


    dEdz = np.zeros_like(E)

    dEdz[1] = -(b_f/v_d)*E[0] + ((a+r_s*b_b)/v_s)*E[1] - ((r_u*b_b)/v_u)*E[2]
    dEdz[2] = (b_b/v_d)*E[0] + ((r_s*b_b)/v_s)*E[1]  - ((a+r_u*b_b)/v_u)*E[2] 

    return dEdz


def ocean_irradiance_analytical(PI, 
                                hbot, 
                                ab_wat, 
                                phy=None, 
                                CDOM_refa=None, 
                                CDOM_chla=None,
                                det=None, 
                                N=30): 
    """
    This is the analytical solution for the three stream irradiance model. The coefficients of 
    absorption are necessarily constant for this method. This means either water only, or a constant profile
    of constituents absorption. 

    """

    ##N centers
    Nm1 = N - 1  

    Ed0 = PI.Ed0
    Es0 = PI.Es0
    Euh = PI.Euh

    coefficients = PI.coefficients
    r_s, r_u, v_d, v_s, v_u = coefficients
 
    ## Default number of shots for BVP shoot method solution
    shots = 3 

    z, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                            ab_wat, 
                                            N, 
                                            Ed0,
                                            coefficients,
                                            phy=phy, 
                                            CDOM_refa=CDOM_refa, 
                                            CDOM_Chla=CDOM_Chla,
                                            det=det, 
                                            grid=PI.grid, 
                                            pt1_perc_zbot=PI.pt1_perc_zbot, 
                                            pt1_perc_phy=PI.pt1_perc_phy)

    ## [For the analytical solution ensure that the coefficients are constant profiles.]
    assert ~np.any(a != a[Nm1])
    assert ~np.any(b != b[Nm1])

    ## [The attenuation depth as calculated from zbot function.]
    Ha = z[0]

    ## [Make the coeffcients into scaalars.]
    a = a[Nm1]
    b = b[Nm1]
    b_b = b_b[Nm1]
    b_f = b_f[Nm1]

    
    ## [Analytical Solution Start.]
    Cs = (a+r_s*b_b)/v_s
    Bu = (r_u*b_b)/v_u
    Fd = b_f/v_d
    Bs = (r_s*b_b)/v_s
    Cu =  (a+r_u*b_b)/v_u
    Bd = b_b/v_d
    Cd = (a+b)/v_d

    D = 0.5*((Cs+Cu) + np.sqrt((Cu+Cs)**2 - 4*Bs*Bu))

    x = -(1/(-(Cs-Cd)*(Cu+Cd) + Bs*Bu)) * ((Cu+Cd)*Fd + Bu*Bd)
    y = -(1/(-(Cs-Cd)*(Cu+Cd) + Bs*Bu)) * (Bs*Fd + Bd*(Cs-Cd))

    ## [The eigenvalues.]
    ## [Lambda minnus. also lambda one in notes.]
    lamm = D - Cu 
    lamp = Cs - D

    ## [Ed at H]
    Edh = analytical_Ed_3stream(Ha, Ed0, a, b, coefficients)

    ## [The numerator of c2.]
    c2num = (Es0*np.exp(lamm*Ha)*(Bs/D) - x*Ed0*np.exp(lamm*Ha)*(Bs/D) + y*Edh)
    c2den = ((Bu/D)*(Bs/D)*np.exp(lamm*Ha) - np.exp(lamp*Ha))
    c2 = c2num / c2den

    c1 = Es0 - x*Ed0 - c2*(Bu/D)

    ## [The solution for Ed]
    Ed = analytical_Ed_3stream(z, Ed0, a, b, coefficients)

    ## [The solution for Es.]
    Es = c1*np.exp(lamm*z) + c2*np.exp(lamp*z)*(Bu/D) + x*Ed

    ## [The solution for Eu.]
    Eu = c1*np.exp(lamm*z)*(Bs/D) + c2*np.exp(lamp*z) + y*Ed

    return Ed, Es, Eu, z


def ocean_irradiance_scipy(PI, 
                           hbot, 
                           ab_wat, 
                           phy = None, 
                           CDOM_refa = None, 
                           CDOM_chla = None, 
                           det = None, 
                           N = 30, 
                           zabb_b=None): 
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
        
        a_r = np.interp(z,zarr,a)
        b_r = np.interp(z,zarr,b)
        b_b_r = np.interp(z,zarr,b_b)
        b_f_r = np.interp(z,zarr,b_f)
        
        dEdz = np.zeros_like(E)

        ## [The independent solution for Ed.]
        dEdz[0] = dEddz_3stream(z, E[0], a_r, b_r, b_b_r, b_f_r, coefficients)
        
        ## [The Es and Eu streams.]
        dEdz[1:] = dEdz_3stream(z, E, a_r, b_r, b_b_r, b_f_r, coefficients)[1:]

        return dEdz
        
    def Ebcs(E_at_h, E_at_0):
        
        return np.array([E_at_0[0] - Ed0, E_at_0[1] - Es0, E_at_h[2] - Euh])

    Ed0 = PI.Ed0
    Es0 = PI.Es0
    Euh = PI.Euh

    coefficients = PI.coefficients

    if zabb_b:  
        zarr, a, b, b_b = zabb_b

        ## [construct the vertical grid to the designated depth.]
        grid = PI.grid
        if grid == 'log': 
            ## log transformed z grid.
            z = Log_Trans_Grid(zarr[0], N) 
        elif grid == 'linear': 
            ## linear z 
            z = np.linspace(zarr[0], 0, N)
        else: 
            raise ValueError("Invalid grid keyword, either 'log' or 'linear'")
 
        a = np.interp(z, zarr, a)
        b = np.interp(z, zarr, b)
        b_b = np.interp(z, zarr, b_b)
        zarr = z
        b_f = b - b_b
    
    else:
        zarr, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                                ab_wat, 
                                                N, 
                                                Ed0,
                                                coefficients,
                                                phy=phy, 
                                                CDOM_refa=CDOM_refa, 
                                                CDOM_chla=CDOM_chla, 
                                                det=det, 
                                                grid=PI.grid, 
                                                pt1_perc_zbot=PI.pt1_perc_zbot, 
                                                pt1_perc_phy=PI.pt1_perc_phy)
 
    Eguess = np.full((3, N), 1.0)

    res = integrate.solve_bvp(derivEdz, Ebcs, zarr, Eguess)

    y = res.y
    z_mesh = res.x

    Ed = np.interp(zarr, z_mesh, y[0])
    Es = np.interp(zarr, z_mesh, y[1])
    Eu = np.interp(zarr, z_mesh, y[2])
    
    return Ed, Es, Eu, zarr


def ocean_irradiance_shootdown(PI, 
                                hbot, 
                                ab_wat, 
                                phy = None, 
                                CDOM_refa = None, 
                                CDOM_chla = None, 
                                det = None,
                                N = 30):
    
    
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
    ##N centers
    Nm1 = N - 1  

    Ed0 = PI.Ed0
    Es0 = PI.Es0
    Euh = PI.Euh

    coefficients = PI.coefficients
 
    ##initial guess... doesn't matter too much
    init_guess = Es0

    E = np.full((N,3), init_guess)
    
    ## BCs included in the initial guess arrays
    E[Nm1,0] = Ed0 ##surface 
    E[Nm1,1] = Es0 ##surface 
    E[0,2] = Euh ##bottom
    
    ## Default number of shots for BVP shoot method solution
    shots = 3 

    z, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                            ab_wat, 
                                            N, 
                                            Ed0,
                                            coefficients,
                                            phy=phy, 
                                            CDOM_refa=CDOM_refa, 
                                            CDOM_chla=CDOM_chla, 
                                            det=det, 
                                            grid=PI.grid, 
                                            pt1_perc_zbot=PI.pt1_perc_zbot, 
                                            pt1_perc_phy=PI.pt1_perc_phy)

    ## [Calculate the independent solution for Ed]
    E[:,0] = numerical_Ed_3stream(z, Ed0, a, b, coefficients)
    
    Eu0_tried = []
    Fmetric = []
     
    for jsh in range(shots) :
    # Integrate down from the top to ensure Ed(1) and Es(1) are good.
 
        if jsh == 0:
            dEu = 0 #-Eu[Nm1]
        elif jsh == 1:
        # for the first case, need some adjustment to get gradient.
            dEu = max(0.01,0.03*E[Nm1,1])
            # dEu = .2
        else: 
            Jslope = (Fmetric[jsh-2]-Fmetric[jsh-1]) / (Eu0_tried[jsh-2]-Eu0_tried[jsh-1]) 
            dEu = -Fmetric[jsh-1]/Jslope
 
#        dEu = max(-E[Nm1,2], dEu)
#        dEu = min(1-E[Nm1,2], dEu)
        E[Nm1,2] = E[Nm1,2] + dEu
        Eu0_tried.append(E[Nm1,2])

        ## [Solving the initial value problem here.]
        for k in range(Nm1, 0, -1): 
            dz = z[k-1] - z[k]
            E[k-1,1:] = RK45(dEdz_3stream, dz, z[k], E[k,:], a[k], b[k], b_b[k], b_f[k], coefficients)[1:]
             
        Fmetric.append(E[0,2] - Euh)

    Ed = E[:,0]
    Es = E[:,1]
    Eu = E[:,2]

    return Ed, Es, Eu, z


def ocean_irradiance_shootup(PI, 
                             hbot, 
                             ab_wat, 
                             phy = None, 
                             CDOM_refa = None, 
                             CDOM_chla = None, 
                             det = None,
                             N = 30):
    
    
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
    ##N centers
    Nm1 = N - 1  

    Ed0 = PI.Ed0
    Es0 = PI.Es0
    Euh = PI.Euh

    coefficients = PI.coefficients
 
    ##initial guess... doesn't matter too much
    init_guess = Es0

    E = np.full((N,3), init_guess)
    
    ## BCs included in the initial guess arrays
    E[Nm1,0] = Ed0 ##surface 
    E[Nm1,1] = Es0 ##surface 
    E[0,2] = Euh ##bottom
    
    ## Default number of shots for BVP shoot method solution
    shots = 3 

    z, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                            ab_wat, 
                                            N, 
                                            Ed0,
                                            coefficients,
                                            phy=phy, 
                                            CDOM_refa=CDOM_refa, 
                                            CDOM_chla=CDOM_chla, 
                                            det=det, 
                                            grid=PI.grid, 
                                            pt1_perc_zbot=PI.pt1_perc_zbot, 
                                            pt1_perc_phy=PI.pt1_perc_phy)

    ## [Calculate the independent solution for Ed]
    E[:,0] = numerical_Ed_3stream(z, Ed0, a, b, coefficients)
    
    Esh_tried = []
    Fmetric = []
     
    for jsh in range(shots) :
    # Integrate down from the top to ensure Ed(1) and Es(1) are good.
 
        if jsh == 0:
            dEs = 0 #-Eu[Nm1]
        elif jsh == 1:
        # for the first case, need some adjustment to get gradient.
            dEs = max(0.01,0.03*E[Nm1,1])
        else: 
            Jslope = (Fmetric[jsh-2]-Fmetric[jsh-1]) / (Esh_tried[jsh-2] - Esh_tried[jsh-1]) 
            dEs = -Fmetric[jsh-1]/Jslope
 
#        dEu = max(-E[Nm1,2], dEu)
#        dEu = min(1-E[Nm1,2], dEu)
        E[0,1] = E[0,1] + dEs
        Esh_tried.append(E[0,1])

        ## [Solving the initial value problem here.]
        for k in range(Nm1): 
            
            dz = z[k+1] - z[k]
            E[k+1,1:] = RK45(dEdz_3stream, dz, z[k], E[k,:], a[k], b[k], b_b[k], b_f[k], coefficients)[1:]
             
        Fmetric.append(E[Nm1,1] - Es0)

    Ed = E[:,0]
    Es = E[:,1]
    Eu = E[:,2]

    return Ed, Es, Eu, z



def ocean_irradiance_semianalytic_inversion(PI, 
                                            hbot, 
                                            ab_wat, 
                                            phy=None, 
                                            CDOM_refa=None, 
                                            CDOM_chla=None, 
                                            det=None, 
                                            N=30):
    
    
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
    coefficients = PI.coefficients
    r_s, r_u, v_d, v_s, v_u = coefficients
    
    Ed0 = PI.Ed0
    Es0 = PI.Es0
    Euh = PI.Euh

    N = N+1
    ##N centers
    Nm1 = N - 1  

    ## [The absorption and scattering calcs.]
    z, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                            ab_wat, 
                                            N, 
                                            Ed0,
                                            coefficients,
                                            phy=phy, 
                                            CDOM_refa=CDOM_refa, 
                                            CDOM_chla=CDOM_chla, 
                                            det=det, 
                                            grid=PI.grid, 
                                            pt1_perc_zbot=PI.pt1_perc_zbot, 
                                            pt1_perc_phy=PI.pt1_perc_phy)
    z = np.flip(z)
    a = np.flip(a)
    b = np.flip(b)
    b_b = np.flip(b_b)
    b_f = np.flip(b_f)
#    z_out = np.zeros(Nm1)
#    for k in range(Nm1):   
#        dz = z[k+1] - z[k]  
#        z_out[k] = z[k] + dz/2 
#    z_out = np.flip(z_out)
 
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
    E_d =  np.flip(numerical_Ed_3stream(np.flip(z), Ed0, np.flip(a), np.flip(b), coefficients))
    
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
    #    ## [The absorption and scattering calcs.]
    z_out, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                            ab_wat, 
                                            N-1, 
                                            Ed0,
                                            coefficients,
                                            phy=phy, 
                                            CDOM_refa=CDOM_refa, 
                                            CDOM_chla=CDOM_chla, 
                                            det=det, 
                                            grid=PI.grid, 
                                            pt1_perc_zbot=PI.pt1_perc_zbot, 
                                            pt1_perc_phy=PI.pt1_perc_phy)
    z_out = np.flip(z_out)
    a = np.flip(a)
    b = np.flip(b)
    b_b = np.flip(b_b)
    b_f = np.flip(b_f)

    Ed = numerical_Ed_3stream(np.flip(z_out), Ed0, np.flip(a), np.flip(b), coefficients)
    Es = np.flip(E_s_z(z_out, z, c_p, c_m, E_d))
    Eu = np.flip(E_u_z(z_out, z, c_p, c_m, E_d))
    z =np.flip(z_out)

    
    #Es = x_lu[:N]
    #Eu = x_lu[N:] 
    #Ed = analytical_Ed(zarr, c)
    
    return Ed, Es, Eu, z


def ocean_irradiance_semianalytic_inversion_ROMS(PI, 
                                            hbot, 
                                            ab_wat, 
                                            phy=None, 
                                            CDOM_refa=None, 
                                            det=None, 
                                            N=30):
    
    """

    This is the Dutkiewicz (2015) Solution found in Appendix B of her paper
    """

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
    coefficients = PI.coefficients
    r_s, r_u, v_d, v_s, v_u = PI.coefficients

    Ed0 = PI.Ed0
    Es0 = PI.Es0
    Euh = PI.Euh
    
    ##N centers
    Nm1 = N - 1  

    ## [The absorption and scattering calcs.]
    zarr, a, b, b_b, b_f = Calc_Abscat_Grid(hbot, 
                                            ab_wat, 
                                            N, 
                                            Ed0,
                                            coefficients,
                                            phy=phy, 
                                            CDOM_refa=CDOM_refa, 
                                            det=det, 
                                            grid=PI.grid, 
                                            pt1_perc_zbot=PI.pt1_perc_zbot, 
                                            pt1_perc_phy=PI.pt1_perc_phy)
 
    c_d = (a+b)/v_d ##used for analytical 
    
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
    
    
    #z = np.linspace(zarr[Nm1], 0, Nm1) ##z array for E_s_z and E_u_z 
    z = zarr
    
    Ed_zarr = numerical_Ed(zarr, c_d, Ed0)
    Es = E_s_z(Nm1, Ed_zarr, z, zarr, c_p, c_m)
    Eu = E_u_z(Nm1, Ed_zarr, z, zarr, c_p, c_m)
    Ed = numerical_Ed(z, c_d, Ed0)
    
    
    return Ed, Es, Eu, z



def ocean_irradiance(PI, 
                     hbot, 
                     ab_wat, 
                     method='scipy', 
                     phy=None,
                     CDOM_refa=None, 
                     CDOM_chla=None, 
                     det=None,
                     N=30, 
                     zabb_b=None): 
    """
    This function runs the ocean irradiance algorithm for the given method.
    """

    if method=='shootdown': 
        Ed, Es, Eu, z = ocean_irradiance_shootdown(PI, hbot, ab_wat, phy=phy, CDOM_refa=CDOM_refa, CDOM_chla=CDOM_chla, det=det, N=N)
    elif method=='shootup': 
        Ed, Es, Eu, z = ocean_irradiance_shootup(PI, hbot, ab_wat, phy=phy, CDOM_refa=CDOM_refa, CDOM_chla=CDOM_chla, det=det, N=N)
    elif method=='scipy':
        Ed, Es, Eu, z = ocean_irradiance_scipy(PI, hbot, ab_wat, phy=phy, CDOM_refa=CDOM_refa, CDOM_chla=CDOM_chla, det=det, N=N, zabb_b = zabb_b)
    elif method=='dutkiewicz':
        Ed, Es, Eu, z = ocean_irradiance_semianalytic_inversion(PI, hbot, ab_wat, phy=phy, CDOM_refa=CDOM_refa, CDOM_chla=CDOM_chla, det=det, N=N)
    elif method=='analytical': 
        Ed, Es, Eu, z = ocean_irradiance_analytical(PI, hbot, ab_wat, phy=phy, CDOM_refa=CDOM_refa, CDOM_chla=CDOM_chla, det=det, N=N)
    else: 
        raise ValueError('Ocean irradiance method unrecognized. Please use one of the follwing methods: \n shootup \n shootdown \n scipy \n dutkiewicz \n analytical')
        

    return Ed, Es, Eu, z
    

def artificial_phy_prof(z,loc,width,conc, prof_type='gauss'):

    if prof_type == 'tan': 
        prof= conc*(1 + np.tanh((z-loc)/width)) 
    elif prof_type == 'gauss': 
        prof = conc* np.exp(-((z-loc)**2)/(2*width**2))
     
    return prof


def Demo(): 
    import matplotlib.pyplot as plt
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    PI = Param_Init()

    PI.pt1_perc_phy = True
    PI.pt1_perc_zbot = True

    PI.grid = 'log'
    
    N = 60
    Nm1 = N-1 
    wavelengths = [443, 551]

    hbot = -700
    
    z_phy = np.linspace(hbot,0,N)

    phy_type = 'Syn'

    phy_prof = artificial_phy_prof(z_phy, -10, 20, 10, prof_type = 'gauss')
    phy_prof = np.full(N, 1)

    fig, axes = plt.subplots(1, 3, sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    
    for k, lam in enumerate(wavelengths):

        ab_wat = abscat(lam, 'water')
        
        ## Define the Phytoplankton class.
        phy = Phy(z_phy, phy_prof, esd(phy_type), abscat(lam, phy_type)[0], abscat(lam, phy_type)[1])

#        Ed, Es, Eu, z = ocean_irradiance_scipy(PI, hbot, ab_wat, phy = phy, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_shootup(PI, hbot, ab_wat, phy=phy, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_shootdown(PI, hbot, ab_wat, phy=phy, N=N)
        Ed, Es, Eu, z = ocean_irradiance_semianalytic_inversion(PI, hbot, ab_wat, phy=phy, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_analytical(PI, hbot, ab_wat, phy=phy, N=N)
 

        markers = ['-.', '-.'] 
        Ed_c = 'g'
        Es_c = 'b'
        Eu_c = 'r'
        if lam == 443:
            ax2.plot(Ed, z, label=f'Ed', color = Ed_c, ls=markers[k] )
            ax2.plot(Es, z, label=f'Es', color = Es_c, ls=markers[k])
            ax2.plot(Eu, z, label=f'Eu', color = Eu_c, ls= markers[k])
        if lam == 551:
            ax3.plot(Ed, z, label=f'Ed', color = Ed_c, ls=markers[k] )
            ax3.plot(Es, z, label=f'Es', color = Es_c, ls=markers[k])
            ax3.plot(Eu, z, label=f'Eu', color = Eu_c, ls= markers[k])

    ax1.plot(phy_prof, z_phy)

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
    
    
    return Ed, Es, Eu, z



def Demo2(): 
    import matplotlib.pyplot as plt
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat
    from ocean_irradiance_module.absorbtion_and_scattering_coefficients import equivalent_spherical_diameter as esd
    from ocean_irradiance_module.PARAMS import Param_Init 
    
    PI = Param_Init()

    PI.pt1_perc_phy = True
    PI.pt1_perc_zbot = True

    PI.grid = 'log'

    PI.Ed0 = 0.7
    PI.Es0 = 0.3
    
    N = 30
    Nm1 = N-1 
    wavelengths = PI.wavelengths

    hbot = -100
    
    z_phy = np.linspace(hbot,0,N)

    phy_type = 'Syn'

    phy_prof = artificial_phy_prof(z_phy, -10, 20, 5, prof_type = 'gauss')
#    phy_prof = np.full(N, 1)

    fig, axes = plt.subplots(1, 5, sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]
    ax5 = axes[4]
    

    for j, lam in enumerate(wavelengths):
        ab_wat = abscat(lam, 'water')
            
        ## Define the Phytoplankton class.
        phy = Phy(z_phy, phy_prof, esd(phy_type), abscat(lam, phy_type)[0], abscat(lam, phy_type)[1])
        
#        cdom_chla = CDOM_chla(z_phy, phy_prof, lam)
        cdom_chla = None

        Ed, Es, Eu, z = ocean_irradiance_scipy(PI, hbot, ab_wat, phy = phy, CDOM_chla=cdom_chla, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_shootup(PI, hbot, ab_wat, phy=phy, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_shootdown(PI, hbot, ab_wat, phy=phy, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_semianalytic_inversion(PI, hbot, ab_wat, phy=phy, N=N)
#        Ed, Es, Eu, z = ocean_irradiance_analytical(PI, hbot, ab_wat, phy=phy, N=N)

        print(Eu[-1])
 
        rgb = W2RGB.wavelength_to_rgb(lam)

#        ax2.plot(Ed, z, color = rgb, ls='-', marker='o', fillstyle='none', label=lam, markersize =3)
#        ax3.plot(Es, z, color = rgb, ls='-', marker='o', fillstyle='none', markersize =3)
#        ax4.plot(Ed +Es, z, color = rgb, ls='-', marker='o', fillstyle='none', markersize =3)
#        ax5.plot(Eu, z, color = rgb, ls='-', marker='o', fillstyle='none', markersize =3)


        ax2.plot(Ed, z, color = rgb, label=lam, ls='-')
        ax3.plot(Es, z, color = rgb, ls='-')
        ax4.plot(Ed +Es, z, color = rgb, ls='-')
        ax5.plot(Eu, z, color = rgb, ls='-')
    ax1.plot(phy_prof, z_phy)

    ax1.set_xlabel('Concentration [mg Chl-a m^-3]')
    ax1.set_ylabel('Z [m]')
    ax1.set_title(r'a) $\rho_{\mathrm{phy}}$')
    ax1.grid()
 
#    ax2.set_xlim(-0.1,1)
    ax2.set_xlabel('Normalized Irradiance')
    ax2.set_title(r'b) $E_d$')
    ax2.legend(title='Wavelength [nm]')
    ax2.grid()

#    ax3.set_xlim(-0.1,1)
    ax3.set_xlabel('Normalized Irradiance')
    ax3.set_title(r'c) $E_s$')
    ax3.grid()

    ax4.set_xlabel('Normalized Irradiance')
    ax4.set_title(r'd) $E_d + E_s$')
    ax4.grid()

#    ax4.set_xlim(-0.1,1)
    ax5.set_xlabel('Normalized Irradiance')
    ax5.set_title(r'e) $E_u$')
    ax5.grid()
    
    fig.show()
    
    
    return Ed, Es, Eu, z


#-------------------------------------MAIN-------------------------------------

if __name__ == '__main__':
    
    Ed, Es, Eu, z = Demo2()
    

