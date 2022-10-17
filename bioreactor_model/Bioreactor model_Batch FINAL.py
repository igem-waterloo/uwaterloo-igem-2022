"""
Waterloo iGEM 2022
Batch Bioreactor Simulation
"""

# Importing Libraries
from matplotlib import colors 
from matplotlib.ticker import PercentFormatter 
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.integrate import solve_ivp
import scipy.optimize as opt

# Creates the growth curve for batch bioreactor

#Design parameters
mu_max = 0.41 #hours^-1, experimental 
Ks = 99 #mg/L, literature
Y_xs = 0.5 #theoretical
Y_ps = 0.15 #approximated
a = 0.4 # growth association parameter, approx.
time = 12 #hours^-1

#Initial conditions
Ao = 90.304 #biomass (X), mg
Bo = 20000 #substrate (S), mg
Co = 0 #product (P), mg

def mu(B):
    ''' Calculates growth rates'''
    growth_rate = mu_max*B/(Ks + B)
    return growth_rate
def monod_eq(mu, Ks, Y_xs, Y_ps, a, time):

    def f(t,y):
        '''ODES'''
        A = y[0] #biomass (X)
        B = y[1] #substrate (S)
        C = y[2] #product (P)
        
        dA_dt = A*mu_max*B/(Ks + B)
        dB_dt = -dA_dt/Y_xs - a*dA_dt/Y_ps
        dC_dt = a*dA_dt

        return np.array([dA_dt, dB_dt, dC_dt]) 

    '''Solves ODES'''
    t_span = np.array([0,time])
    times = np.linspace(t_span[0], t_span[1], time)
    #times = ([0,1.92,3.5,9.25,13.75])
    y0 = np.array([Ao,Bo,Co])
    soln = solve_ivp(f, t_span, y0, t_eval = times)
    
    '''Creates Dataframe'''
    t = soln.t
    A = soln.y[0]
    B = soln.y[1]
    C = soln.y[2]
    dataa = [t,A,B,C]
    dataT = np.transpose(dataa)
    df = pd.DataFrame(data = dataT)
    df.columns = ['Time (hours)','Biomass (X)', 'Substrate (S)', 'Product (P)']
 
    print(df)
    return df

def testing(mu, Ks, Y_xs, Y_ps, a, time):
    x_time = ([0,1.92,3.5,9.25,13.75])
    y_real = ([90.304, 118.192, 233.949, 1207.816, 1424.722]) #biomass exp.
    
    testing = monod_eq(mu, Ks, Y_xs, Y_ps, a, time)
    #plt.legend(['Time','Biomass', 'Substrate', 'Product'])
    plt.plot(testing['Biomass (X)'])
    plt.plot(x_time, y_real)
    plt.legend(['Theoretical', 'Real'])
    plt.xlabel('Time (h)')
    plt.ylabel('Biomass (mg)')
    plt.title('Theroetical vs Real Biomass Production')
    return None
    
testing(mu, Ks, Y_xs, Y_ps, a, time) 
#beginning of muting

def squared_error(par):
    print(par)
    um = par[0]
    Ks = par[1]
    Y_xs = par[2]
    Y_ps = par[3]
    a = par[4]
    
    
    training_set = pd.read_csv('BioreactorTest.csv')
    #df2 = pd.DataFrame(data = [cell_dens]*10)
    test_set = monod_eq(par[0], par[1], par[2], par[3], par[4],time)
    
    
    sse = 0
    for i,row in training_set.iterrows():
        sse = sse + (training_set.iloc[(i,0)] - test_set.iloc[(i,0)])**2
    print(sse)
        
    return sse
     
   # Optimization Funciton 

#training_set = monod_eq(mu, Ks, Y_xs, Y_ps, a, time)

par_all = [mu_max, Ks, Y_xs, Y_ps, a]


# Test points are defined as monod equaion
testpoints = monod_eq(mu, Ks, Y_xs, Y_ps, a, time)
s_error = lambda p: squared_error(p)
solution = opt.minimize(s_error, par_all, method = 'Nelder-Mead', tol = 1e-6)

plt.figure()
plt.plot(testpoints)
plt.title('Batch Simulation')
plt.xlabel("Time (h)")
plt.ylabel('Mass of X,S,P (mg)')
plt.legend(['Time','Biomass', 'Substrate', 'Product'])

