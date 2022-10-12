"""
Waterloo iGEM 2022
Fed Batch Bioreactor Simulation
"""

# Importing Libraries
from matplotlib import colors 
from matplotlib.ticker import PercentFormatter 
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.integrate import solve_ivp
import scipy.optimize as opt

#Defining parameters
mu_max = 0.4 #hours^-1, experimental 
Ks = 99 #mg/L, literature val
Y_xs = 0.5 #experimental 
Y_ps = 0.15 #approximation
a = 0.4 #growth association const., approximation
time = 20 #hours^-1
Feed = 1 # input flowrate, L/h
Bf = 40000 # input concentration, mg/L

#Initial conditions
Ao = 90.304 #mg, biomass (X) 
Bo = 20000 #mg, substrate (S) 
Co = 0 #mg, product (P)
Vol = 5 # Initial volume, L 

def mu(B):
    ''' Calculates growth rates'''
    growth_rate = mu_max*B/(Ks + B)
    return growth_rate

def F(t):
    '''Feed flowrate'''
    return Feed

def monod_eq(mu, Ks, Y_xs, Y_ps, a, time):
    '''Calculates growth curves'''
    
    def f(t,y):
        '''ODES'''
        A = y[0] #biomass (X)
        B = y[1] #substrate (S)
        C = y[2] #product (P) 
        V = y[3] #volume of reactor
    
        dV_dt = F(t)
        dA_dt =(A*mu(B))-(dV_dt*A/V)
        dB_dt = (F(t)*(Bf-B)/V) - (A*mu(B)/Y_xs) - (a*A*mu(B)/Y_ps)
        dC_dt = a*A*mu(B)-(dV_dt*C/V)

        return np.array([dA_dt, dB_dt, dC_dt, dV_dt]) 

    '''Solves ODES'''
    t_span = np.array([0,time])
    times = np.linspace(t_span[0], t_span[1], time)
    y0 = np.array([Ao,Bo,Co,Vol])
    soln = solve_ivp(f, t_span, y0, t_eval = times)
    
    '''Creates Dataframe'''
    t = soln.t
    A = soln.y[0]
    B = soln.y[1]
    C = soln.y[2]
    V = soln.y[3]
    dataa = [A,B,C,V]
    dataT = np.transpose(dataa)
    df = pd.DataFrame(data = dataT)
    df.columns = ['Biomass (X)', 'Substrate (S)', 'Product (P)', 'Volume (V)']
    
    print(df)
    return df

testpoints = monod_eq(mu, Ks, Y_xs, Y_ps, a, time)
testpoints2 = testpoints.loc[:,['Biomass (X)', 'Substrate (S)', 'Product (P)']]

corrected_testpoints = testpoints
for i in range(1,time):
    if testpoints.loc[i, 'Substrate (S)'] < 0:
        corrected_testpoints.loc[i, 'Substrate (S)'] = 0
        corrected_testpoints.loc[i, 'Biomass (X)'] = testpoints.loc[i-1, 'Biomass (X)'] #+ testpoints.loc[i-1, 'Biomass (X)']*mu(3)
        corrected_testpoints.loc[i, 'Product (P)'] = testpoints.loc[i-1, 'Product (P)'] #+ testpoints.loc

corrected_testpoints2 = corrected_testpoints.loc[:,['Biomass (X)', 'Substrate (S)', 'Product (P)']]
lastrow = time - 1
P_final = round(testpoints.loc[lastrow,'Product (P)'],2)

plt.plot(corrected_testpoints2)
plt.xlabel("Time (h)")
plt.ylabel('Mass of X,S,P (mg)')
plt.legend(['Biomass', 'Substrate', 'Product'], loc = "upper right")
plt.title(f'Feed Concentration = {Bf}mg/L ; V = {Vol}L')
plt.suptitle("Fed-Batch Simulation")
#plt.text(13,18000,f" P final = {P_final}mg")
plt.text(13,21000,f" P final = {P_final}mg",bbox=dict(facecolor='white', alpha=0.2))
