
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
mu_max = 0.4 #hours^-1, experimental (found in textbook, pg 187)
Ks = 99 #mg/L, experimental (found in textbook, pg 187)
Y_xs = 0.5 #experimental 
Y_ps = 0.15 #need lit value
a = 0.4 # need lit value
time = 40#hours^-1
F = 1 # input flowrate, L/h
Bf = 6000# mg/L concentration pure substate in as feed 
V = 3#liter 


#Initial conditions
Ao = 90.304 #biomass (X) (should be mg?)
Bo = 20000 #substrate (S) (should be mg?)
Co = 0 #product (P) (should be mg?)
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
    times = np.linspace(t_span[0], t_span[1],time)
    y0 = np.array([Ao,Bo,Co,Vol])
    
    soln = solve_ivp(f, t_span, y0, t_eval = times)
    t = soln.t
    A = soln.y[0]
    B = soln.y[1]
    C = soln.y[2]
    #V = soln.y[3]
 
    dataa = [A,B,C]
    dataT = np.transpose(dataa)

    df = pd.DataFrame(data = dataT)
    df.columns = ['Biomass (X)', 'Substrate (S)', 'Product (P)']
    #df['B'] = B
    #df['C'] = C
    print(df)
    return df


# Test points are defined as monod equaion
testpoints = monod_eq(mu, Ks, Y_xs, Y_ps, a, time)

corrected_testpoints = testpoints
for i in range(1,time):
    if testpoints.loc[i, 'Substrate (S)'] < 0:
        corrected_testpoints.loc[i, 'Substrate (S)'] = 0
        corrected_testpoints.loc[i, 'Biomass (X)'] = testpoints.loc[i-1, 'Biomass (X)'] #+ testpoints.loc[i-1, 'Biomass (X)']*mu(3)
        corrected_testpoints.loc[i, 'Product (P)'] = testpoints.loc[i-1, 'Product (P)'] #+ testpoints.loc

lastrow = time - 1
P_final = round(testpoints.loc[lastrow,'Product (P)'],1)

plt.plot(corrected_testpoints) 
plt.title(f' Feed Concentration = {Bf}mg/L ; V = {V}L')
plt.suptitle("CSTR Simulation")
plt.xlabel("Time (h)")
plt.ylabel('Mass of X,S,P (mg)')
#plt.legend(['Biomass', 'Substrate', 'Product','Volume'])
plt.legend(['Biomass', 'Substrate', 'Product',], loc = "upper right")
plt.text(26.7,27000,f" P final = {P_final}mg",bbox=dict(facecolor='white', alpha=0.2))