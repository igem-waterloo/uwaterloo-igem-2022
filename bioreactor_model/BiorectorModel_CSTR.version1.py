
# Importing Libraries
from matplotlib import colors 
from matplotlib.ticker import PercentFormatter 
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.integrate import solve_ivp
import scipy.optimize as opt

# Creates the growth curve for CSTR  bioractor
#Design parameters
mu_max = 0.4 #hours^-1, experimental 
Ks = 99 #mg/L, literature
Y_xs = 0.5 #experimental 
Y_ps = 0.15 #approx.
a = 0.4 # approx.
time = 15  #hours^-1
F = 1 # input flowrate, L/h
Bf = 20000 # mg/L concentration pure substate in as feed 
V = 4 #tank volume, liter 


#Initial conditions
Ao = 90.304 #biomass (X), mg
Bo = 20000 #substrate (S), mg
Co = 0 #product (P), mg
Vol = V

def mu(B):
    ''' Calculates growth rates'''
    growth_rate = mu_max*B/(Ks + B)
    return growth_rate

def monod_eq(mu, Ks, Y_xs, Y_ps, a, time):
    '''Calculates growth curves'''
    
    def f(t,y):
        '''ODES'''
        A = y[0] #biomass (X)
        B = y[1] #substrate (S)
        C = y[2] #product (P) 
        V = y[3]
        

        dadt = mu(B)*A
		
        # dV_dt = 0 , volume is constant in a CSTR
        dV_dt = 0
        dA_dt = dadt*V - F*A
        dB_dt = F*Bf- F*B - (dadt/Y_xs) - (a*dadt/Y_ps)
		     # F*(Bf) is the flow in of substrate
        dC_dt = a*dadt*V - F*C

        return np.array([dA_dt, dB_dt, dC_dt,dV_dt])

    t_span = np.array([0,time])
    times = np.linspace(t_span[0], t_span[1], time)

    y0 = np.array([Ao,Bo,Co,Vol])
    
    soln = solve_ivp(f, t_span, y0, t_eval = times)
    t = soln.t
    A = soln.y[0]
    B = soln.y[1]
    C = soln.y[2]
    V = soln.y[3]
 
    dataa = [A,B,C,V]
    dataT = np.transpose(dataa)

    df = pd.DataFrame(data = dataT)
    #df['B'] = B
    #df['C'] = C
    print(df)
    return df


# Test points are defined as monod equaion
testpoints = monod_eq(mu, Ks, Y_xs, Y_ps, a, time)

plt.plot(testpoints) 
plt.title(f'CSTR Simulation . Bf={Bf} V={V}')
plt.xlabel("Time (h)")
plt.ylabel('Mass of X,S,P (mg)')
plt.legend(['Biomass', 'Substrate', 'Product','Volume'])