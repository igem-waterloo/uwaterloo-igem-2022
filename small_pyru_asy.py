import numpy as np
from math import log
from scipy import integrate

# the original function has a ~ in place of _, I assume this 
# means that it's not used in the function
def asy_reaction_ODE(_, y, w):
    dy = np.zeros([2,1], dtype=float)
    # again, ./ in matlab is element-wise right division
    # in numpy this is simply /
    dy[1] = -y[1]/(w + y[1]*y[2])
    dy[2] = y[1]/(w + y[1]*y[2]) - y[2]
    return dy

""" 
original lines:
    options = odeset('RelTol',2.3e-14,'AbsTol',2.3e-35)
    [T,Y] = ode15s(@(t,y) asy_reaction_ODE(t, y, w),[Tmin Tmax],y0,options)
create an ODE model with:
    1. relative error tolerance (RelTol in matlab) = 2.3e-14
    2. absolute error tolerance (AbsTol in matlab) = 2.3e-35
    3. 15th order
    4. where asy_reaction_ODE is the initial system
    5. 15s stands for stiff solver
suggested to used this for python
scipy.integrate.ode(initial system).set_integrator('vode', method='bdf', order=15)
    1. bdf: backward differentiation formula -- it's recced to use LSODA which 
        seems to be diff in that it does auto stiffness switching
    2. otherwise works
    3. vode: Real-valued Variable-coefficient Ordinary Differential Equation solver
        seems to be what provides the bdf option
"""
def small_pyru_asy(alp, K_1M, K_1i, l, ep):
    Tmin = -1e4
    Tmax = 50
    w = alp*K_1M
    y0 = [-Tmin - w*log(-Tmin)/2, 1]
    T, Y = integrate.ode(asy_reaction_ODE(t, y, w)).set_integrator('vode', method='bdf', order=15)
    S1 = Y[:1]
    S2 = Y[:2]
    return T, S1, S2