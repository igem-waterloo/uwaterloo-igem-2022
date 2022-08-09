import numpy as np
from math import sqrt, exp, log
from scipy import integrate

import small_pyru_asy

def asymptotic_solutions(T, k, ep):

    k1 = k[1]
    k2 = k[2]*ep
    k_2 = k[3]*ep
    k3 = k[4]
    k4 = k[5]
    k5 = k[6]
    k6 = k[7]
    k7 = k[8]

    K_1M = k[9]
    K_1i = k[10]/ep
    K_2M = k[11]
    K__2M = k[12]/ep
    K_3aM = k[13]/ep # unused
    K_3bM = k[14]/ep
    K_3i = k[15]/ep # unused
    K_4M = k[16]/ep
    K_5M = k[17]
    K_6M = k[18]
    K_7M = k[19]
    A = k[20]

    v2 = k2/K_2M
    v_2 = k_2/K__2M
    v3 = k3/K_3bM
    v4 = k4/K_4M
    v5 = k5/K_5M
    v6 = k6/K_6M
    v7 = k7/K_7M

    alp = A + v2*v3/v_2

    y_key = (v2*v3/v_2)*sqrt(K_1i*(1 - exp(-2*alp*T))/alp)

    S1 = 1 - sqrt(ep*alp*K_1i)*(T + log(1 + sqrt(1 - exp(-2*alp*T)))/alp)
    S2 = ep**(1/2)*(v_2/(v2*v3))*y_key
    S3 = ep**(3/2)*y_key/v3
    S4 = ep**(3/2)*y_key/v4

    S5 = np.zeros(len(S1), dtype=float)
    S6 = S5
    S7 = S5

    # this is to select elements in list which are less than 10, then greater than 10
    lower = [i for i in range(len(T)) if T[i] < 10]
    higher = [i for i in range(len(T)) if T[i] > 10]

    t_int = T[lower]

    # * in matlab is element-wise multiplication, it seems that the python
    # equivalent is just * for the whoole matrix
    integral5 = ep**(1/2)*(v2*v3/(v_2))*sqrt(K_1i/alp)*exp(-v5*t_int) \
        *integrate.cumtrapz(t_int,exp(v5*t_int)*sqrt(1 - exp(-2*alp*t_int)))
    S5[lower] = integral5
    S5[higher] = integral5[-1]

    integral6 = exp(-v6*t_int) \
        *integrate.cumtrapz(t_int,exp(v6*t_int)*v5*S5[lower])
    S6[lower] = integral6
    S6[higher] = integral6[-1]

    integral7 = exp(-v7*t_int) \
        *integrate.cumtrapz(t_int,exp(v7*t_int)*v6*S6[lower])
    S7[lower] = integral7
    S7[higher] = integral7[-1]

    P = integrate.cumtrapz(T, v7*S7)

    S = [S1,S2,S3,S4,S5,S6,S7,P]

    # changed lambda to l, since it's a keyword in python
    l = K_1i*(v2*v3*(1/K_2M - v2/k_2 + v2*v3/(v_2*k3))/(alp*v_2) + alp)/2

    S1_med = 1 - sqrt(alp*K_1i*ep)*T + \
        ep*l*T - sqrt(ep)*sqrt(K_1i/alp)*log(2) - \
        (K_1M/2)*sqrt(alp*K_1i*ep)*(log(1 - sqrt(alp*K_1i*ep)*T))

    [T_long,S1_long_ND,S2_long_ND] \
        = small_pyru_asy(alp,K_1M,K_1i,l,ep)

    t_long = 1/sqrt(ep*alp*K_1i) + T_long/alp

    S1_long = sqrt(ep*K_1i/alp)*S1_long_ND
    S2_long = sqrt(ep*K_1i/alp)*S2_long_ND
    S3_long = ep*S2_long*v2/v_2
    S4_long = S3_long*v3/v4
    S5_long = ep**(1/2)*(v2*v3/(v_2))*sqrt(K_1i/alp)*exp(-v5*T_long) \
        *integrate.cumtrapz(T_long,exp(v5*T_long)*S2_long_ND)
    S6_long = v5*exp(-v6*T_long) \
        *integrate.cumtrapz(T_long,exp(v6*T_long)*S5_long)
    S7_long = v6*exp(-v7*T_long) \
        *integrate.cumtrapz(T_long,exp(v7*T_long)*S6_long)
    P_long = integrate.cumtrapz(T_long,v7*S7_long)
    S_long = [S1_long,S2_long,S3_long,S4_long,S5_long,S6_long,S7_long,P_long]

    return S, t_long, S_long, S1_med