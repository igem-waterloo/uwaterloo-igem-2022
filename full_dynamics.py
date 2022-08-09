# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 23:58:13 2022

@author: rxche
"""

import parameters as p
import scipy
import reaction_ODE as reaction_ode
import time

def full_dynamics():
    """
    Main function to solve the system of ODE's for determining
    IDP concentration
    Returns
    -------
    None.

    """
    
    dim = "non-dim"
    pyruvate = "infinite"
    
    # Import parameters file
    k, y0, ep = p.params()
    
    # k: reaction rate constant array
    # y0: metabollites array
    # ep: epsilon for the artificial small dimensionsless parameter
    
    s_0 = y0[0] # initial level of pyruvate in system
    
    if dim == "non-dim":
        # Non-dimensionalize the system variables by scaling each metabolite
        # concentration with s_0
        y0 = y0 / s_0
        k[0:7] = k[0:7] / k[0]
        k[8:18] = k[8:18] / s_0
    # If type is other than dimension, raise unexpected type
    elif dim != "dim":
        raise ValueError('Unexpected type. Cannot proceed')
        
    tmax = 30 # temperature of grown bacteria cultures
    # Start timer
    t_start = time.time()
    # Integrate using stiff ODE solver
    t, y = scipy.integrate.ode(reaction_ode(y0, k, s_0, ep, dim, pyruvate)).set_integrator('vode', method='bdf', order=15)
        
    elapsed = time.time() - t_start
    print(f"Elapsed time: {elapsed}")
        
full_dynamics()