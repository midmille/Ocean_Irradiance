"""

Author: Miles D. Miller UCSC 
Created: 12:32pm April 10, 2022 

This file contains ODE schemes.
"""

import numpy as np

def Euler_Forward(N, dt, t, dydt, y0): 
    """
    This method implements the Euler forward scheme for the solution to the given 
    vector valued ODE. This is originally used for AM 213B Homework 1, Question 2. 

    The algorith assumes a uniform grid. 

    It also assumes the y is a vector.

    Parameters
    ----------
    N: Integer
        The number of grid points. 
    dt: Float
        The time step. 
    t: 1-D Array, [N]
        The uniform time grid. 
    dydt: Function
        The rhs side of the differential equation. 
    y0: 1-D Array, [M]
        The initial condition of the IVP. 

    Returns 
    -------
    y: 2-D Array, [N, M]
        The solution to the ODE, with y as a vector. 
 
    """
    ## [The length of y0 tells us how many vectors there are]
    M = len(y0)

    ## [Init the solution y array.]
    y = np.zeros((N,M))

    ## [Set the initial conditions.]
    y[0,:] = y0

    ## [Loop over the t-grid.]
    for k in range(N-1): 
        y[k+1,:] = y[k,:] + dt * dydt(y[k,:], t[k])

    return y


def RK2(N, dt, t, dydt, y0):
    """
    This function implements an explicit 2 stage Runge Kutta method for the solution
    of the given ODE. This is used for the solution to 213B Homework 1, Question 2. 

    The algorithm assumes a uniform grid. 

    It also assumes that y is a vector. 

    Parameters
    ----------
    N: Integer
        The number of grid points. 
    dt: Float
        The time step. 
    t: 1-D Array, [N]
        The uniform time grid. 
    dydt: Function
        The rhs side of the differential equation. 
    y0: 1-D Array, [M]
        The initial condition of the IVP. 

    Returns 
    -------
    y: 2-D Array, [N, M]
        The solution to the ODE, with y as a vector. 
        
    """

    ## [The length of y0 tells us how many vectors there are]
    M = len(y0)

    ## [Init the solution y array.]
    y = np.zeros((N,M))

    ## [Set the initial conditions.]
    y[0,:] = y0

    ## [Loop over the t-grid.]
    for k in range(N-1): 
        
        ## [The three stages of RK3]
        k1 = dydt(y[k,:], t[k])
        k2 = dydt(y[k,:] + dt*k1, t[k] + dt)

        ## [Solving for the next step in y.]
        y[k+1, :] = y[k,:] + 0.5*dt*(k1 + k2)

    return y
   

def RK3(N, dt, t, dydt, y0): 
    """
    This function implements an explicit 3 stage Runge Kutta method for the solution
    of the given ODE. This is used for the solution to Homework 1, Question 1, Part B, Section 2
    in 213B. 

    The algorithm assumes a uniform grid. 

    It also assumes that y is a vector. 

    Parameters
    ----------
    N: Integer
        The number of grid points. 
    dt: Float
        The time step. 
    t: 1-D Array, [N]
        The uniform time grid. 
    dydt: Function
        The rhs side of the differential equation. 
    y0: 1-D Array, [M]
        The initial condition of the IVP. 

    Returns 
    -------
    y: 2-D Array, [N, M]
        The solution to the ODE, with y as a vector. 
        
    """

    ## [The length of y0 tells us how many vectors there are]
    M = len(y0)

    ## [Init the solution y array.]
    y = np.zeros((N,M))

    ## [Set the initial conditions.]
    y[0,:] = y0

    ## [Loop over the t-grid.]
    for k in range(N-1): 
        
        ## [The three stages of RK3]
        k1 = dydt(y[k,:], t[k])
        k2 = dydt(y[k,:] + dt*0.5*k1, t[k] + 0.5*dt)
        k3 = dydt(y[k,:] + dt*(-1*k1 + 2*k2), t[k] + dt)

        ## [Solving for the next step in y.]
        y[k+1, :] = y[k,:] + dt*((1/6)*k1 + (2/3)*k2 + (1/6)*k3)

    return y


def AB3(N, dt, t, dydt, y0):
    """ 
    The Adams Bashforth three step method. Using RK3 for the first three steps.

    The algorithm assumes a uniform grid. 

    It also assumes that y is a vector. 

    Parameters
    ----------
    N: Integer
        The number of grid points. 
    dt: Float
        The time step. 
    t: 1-D Array, [N]
        The uniform time grid. 
    dydt: Function
        The rhs side of the differential equation. 
    y0: 1-D Array, [M]
        The initial condition of the IVP. 

    Returns 
    -------
    y: 2-D Array, [N, M]
        The solution to the ODE, with y as a vector. 
        
    """
    ## [The length of y0 tells us how many vectors there are]
    M = len(y0)

    ## [Init the solution y array.]
    y = np.zeros((N,M))

    ## [Set the initial conditions.]
    y[0,:] = y0

    ## [Use RK3 to initialize the first 3 steps.]
    y[:3,:] = RK3(3, dt, t, dydt, y0)

    ## [Loop over the rest.]
    for k in range(N-3): 
        ## [AB3] 
        y[k+3,:] = y[k+2,:] + dt/12 * (23*dydt(y[k+2,:], t[k+2]) - 16*dydt(y[k+1,:], t[k+1]) + 5*dydt(y[k+1,:], t[k+1]))

    return y


def Midpoint_Implicit(N, dt, t, dydt, y0): 
    """
    This method uses the implicit midpoint method
    The algorithm assumes a uniform grid. 

    It also assumes that y is a vector. 

    Parameters
    ----------
    N: Integer
        The number of grid points. 
    dt: Float
        The time step. 
    t: 1-D Array, [N]
        The uniform time grid. 
    dydt: Function
        The rhs side of the differential equation. 
    y0: 1-D Array, [M]
        The initial condition of the IVP. 

    Returns 
    -------
    y: 2-D Array, [N, M]
        The solution to the ODE, with y as a vector. 
        
    """


    return 
    

        
 
