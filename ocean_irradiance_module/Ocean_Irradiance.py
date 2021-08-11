"""
Created on Thu Dec 10 08:52:28 2020

@author: Miles Miller
"""
import numpy as np

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


def zbot_func(E_d_0, a_wat, b_wat, v_d):
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
    c = (a_wat + b_wat) / v_d
    Ed = analytical_Ed(zbots, c, E_d_0)

    for k, Ed in enumerate(Ed) :
        EdoE_d_0 = Ed / E_d_0
        if EdoE_d_0 < .001 :
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


def ocean_irradiance(hbot, Ed0, Es0, Euh, ab_wat, phy = None, N = 30, pt1_perc_zbot = True):
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
    r_s = 1.5 
    r_u = 3.0 
    
    v_d = .9 
    v_s = .83 
    v_u = .4 
    
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
        zbot_pt1perc = zbot_func(Ed0, a_wat, b_wat, v_d)
        ## choosing the smaller zbot and making negative
        zbot = -min(abs(hbot), abs(zbot_pt1perc))
    elif pt1_perc_zbot == False: 
        zbot = hbot 
    ## log transformed z grid.
    z = Log_Trans(zbot, N) 
    
    ## Interpolating a,b vectors from z_phy to z.
    ## Should I create another z_grid that denotes the centers for the a,b below
    a = np.interp(z,z_phy,a)
    b = np.interp(z,z_phy,b)
        
    # if N != len(Ed1) or N !=len(a)+1 or N !=len(b)+1 :
    #     print('lengths of z and Ed must be the same, and a&b should be 1 less.')
        

    b_b = .551*b 
    b_f = b - b_b 
    
    ##coefficient of downward direct irradiance 
    c_d = (a+b)/v_d 
    
    Ed=np.copy(Ed1)
    Es=np.copy(Es1)
    Eu=np.copy(Eu1)

    Eu0_tried = []
    Fmetric = []
    
    for jsh in range(shots) :
   # Integrate down from the top to ensure Ed(1) and Es(1) are good.

        if jsh == 0:
            dEu = 0
        elif jsh == 1:
        # for the first case, need some adjustment to get gradient.
            dEu = max(0.01,0.03*Es[Nm1])
        else: 
            Jslope = (Fmetric[jsh-2]-Fmetric[jsh-1]) / (Eu0_tried[jsh-2]-Eu0_tried[jsh-1]) 
            dEu = -Fmetric[jsh-1]/Jslope

            
        Eu[Nm1] = Eu[Nm1] + dEu
        Eu0_tried.append(Eu[Nm1])

        Edmid = np.zeros(Nm1)
        ## integrate Es down the water column
        ## The range does not actually go to k=-1, only to k=0. 
        ## i.e. does not include stop point of range. 
        for k in range(Nm1-1 , -1, -1) :
         # carry out a RK4 algorithm
    
         # z is +ve upward, opposite to that in Dutkiewicz
         # steping DOWNWARD through water column means dz should be less than 0
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
            
            ## calculate a metric that indicates goodness of our shot.
            ## since Eu(bot) = 0, our metric of fit is just the bottom value for Eu.
        Fmetric.append(Eu[0])


    return Ed, Es, Eu, z


def artificial_phy_prof(z,loc,width):

    prof = .5*(1 + np.tanh((z-loc)/width))
    
    return prof


def Demo(): 
    ## Demonstration of finding the irradiance profiles for an artificial 
    ## phytoplankton concentration profile
    
    ## Absorbtion and scattering coefficients are taken from Dutkiwicz 2015 
    ## for the Diatom species. 
    import matplotlib.pyplot as plt
    
    N = 301
    Nm1 = N-1 
    
    z = np.linspace(-300,0,N)

    phy_prof = artificial_phy_prof(z, -70, 40)
    
    ab_wat = (.01,.004)
    a_phy = .01
    b_phy =  .0039
    
    ## Define the Phytoplankton class.
    phy = Phy(z, phy_prof, a_phy, b_phy)

    Ed0 = .7
    Es0 = 1 - Ed0
    Euh = 0 

    zbot = z[0]

    Ed, Es, Eu, zarr = ocean_irradiance(zbot,Ed0,Es0,Euh,ab_wat,phy, N=N)

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
    ax2.set_xlabel('Irradiance')
    ax2.set_title('Resulting Irradiance Profiles')
    ax2.legend()
    ax2.grid()
    
    plt.show()
    return 


#-------------------------------------MAIN-------------------------------------

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Ocean Irradiance Fucntion')

    parser.add_argument('--demo', action='store_true',
                        help='Run the demo and plot')
    
    args = parser.parse_args()
    
    if args.demo: 
        Demo()





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
