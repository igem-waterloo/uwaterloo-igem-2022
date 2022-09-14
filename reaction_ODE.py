# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 17:26:31 2022

@author: rxche
"""
import numpy as np

def reaction_ODE(t, y, k, s_0, ep, dim, pyruvate):
    """
    Function containing the ODEs to be solved 

    Parameters
    ----------
    y : Array
        Initial concentration conditions
    k : Array
        Reaction rates
    s_0 : float
        Initial value of pyruvate concentration for dimensionless system
    ep : Float
        Small dimensionless parameter epsilon
    dim: String
        Dimension or dimensionless type
    pyruvate : string
        Finite or infinite type

    Returns
    -------
    None.

    """
    
    # Initialize array to hold ODE values
    dy = np.zeros(len(y), 1)
    
    # Define the enzymes
    
    s1 = y[0] # Pyruvate
    s2 = y[1] # Acetyl CoA
    s3 = y[2] # Acetoacetyl CoA
    s4 = y[3] # HMG CoA
    s5 = y[4] # Mevalonate
    s6 = y[5] # Mevalonate-phosphate
    s7 = y[6] # Mevalonate diphosphate
    s8 = y[7] # Isopentyl diphosphate
    
    # Define the reaction rates.
    k1 = k[0] # k1 saccharomyces cerevisiae
    k2 = k[1] # k2 enterococcus faecalis
    k_2 = k[2] # k_2 enterococcus faecalis
    k3 = k[3] # k3 saccharomyces cerevisiae
    k4 = k[4] # k4 enterococcus faecalis
    k5 = k[5] # k5 methanosarcina mazei
    k6 = k[6] # k6 saccharomyces cerevisiae
    k7 = k[7] # k7 saccharomyces cerevisiae
        
    # Define the dimensionless parameters
    k8 = k[8] # k_1^M saccharomyces cerevisiae
    k9 = k[9] # k_1^i saccharomyces cerevisiae
    k10 = k[10] # k_2^M enterococcus faecalis
    k11 = k[11] # k_{-2}^M enterococcus faecalis
    k12 = k[12] # k_{3,a}^M saccharomyces cerevisiae
    k13 = k[13] # k_{3,b}^M saccharomyces cerevisiae
    k14 = k[14] # k_3^i saccharomyces cerevisiae
    k15 = k[15] # k_4^M enterococcus faecalis
    k16 = k[16] # k_5^M methanosarcina mazei
    k17 = k[17] # k_6^M saccharomyces cerevisiae
    k18 = k[18] # k_7^M saccharomyces cerevisiae
    A = k[19] # A
    
    # Define the reactions
    if dim == "dim":
        # Re-define the initial value of pyruvate concentration and dimensionless parameter
        s_0 = 1
        ep = 1
    # If not a dimension or dimensionless system raise ValueError
    elif dim != "non-dim":
        raise ValueError("Unexpected type. Cannot proceed")
    
    # The following correspond to the ODE's (equations 2 in the paper) - dimensional system
    r1 = s1 / (s1 + k8 + s1 * s2 / k9)
    r2 = k2 * s2 / (s2 + k10)
    r_2 = k2 * s3 / (s3 + k10)
    
    if (s2 == 0) and (s3 == 0):
        r3 = 0
    else:
        r3 = (k3 * s2 * s3) / (s2 * s3 + k12 * s3 * (1 + s3 / k14) + k13 * s2)
    
    r4 = k4 * s4 / (s4 + k15)
    r5 = k5 * s5 / (s5 + k16)
    r6 = k6 * s6 / (s6 + k17)
    r7 = k7 * s7 / (s7 + k18)
    rA = A * s2
    
    # Define the reaction kinetics
    
    # Substrate kinetics
    if pyruvate == "finite":
        dy[0] = -r1
    elif pyruvate == "infinite":
        dy[0] = 0
    else:
        raise ValueError("Unexpected type. Cannot proceed")
    
    # The following correspond to the ODE's (equations 3 in the paper) - dimensionless system
    dy[1] = r1 - r2 + r_2 - rA
    dy[2] = r2 - r_2 - r3
    dy[3] = r3 - r4
    dy[4] = r4 - r5
    dy[5] = r5 - r6
    dy[6] = r6 - r7
    dy[7] = r7
    
    return dy
    