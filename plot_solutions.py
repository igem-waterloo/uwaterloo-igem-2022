import numpy as np
from datetime import time
import matplotlib.pyplot as plt

T,Y,k,ep,pyruvate,Tmax = Full_dynamics

SeaGreen_l = [124,205,124]/255
SeaGreen_d = [84,139,84]/255

SteelBlue_l = [99,184,255]/255
SteelBlue_d = [54,100,139]/255

Red_l = [255,64,64]/255
Red_d = [139,35,35]/255

Choc_l = [255,127,36]/255
Choc_d = [210,105,30]/255

Pink_l = [255,110,180]/255
Pink_d = [139,58,98]/255

Col_l = [SeaGreen_lSteelBlue_lRed_lChoc_lPink_l]
Col_d = [SeaGreen_dSteelBlue_dRed_dChoc_dPink_d]

start = time.time
S,t_long,S_long,S1_med = asymptotic_solutions(T,k,ep)
print(time.time-start)

# size(S)

# S1 = S(1,:)
# S2 = S(2,:)
# S3 = S(3,:)
# S4 = S(4,:)

if pyruvate == 'finite':
    fig, ax = plt.subplots()
    for j in range(4):
        ax.plot(T,Y[:,j], color=Col_l[j,:],linewidth=3, label=f'S{j+1}')
        ax.plot(T,S[:,j],'--', color=Col_d[j,:],linewidth=3)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1e-5, Tmax)
    ax.set_ylim(1e-20, 1)
    plt.legend(loc='best')
        
    fig.show()

    fig, ax = plt.subplots()
    for j in range(4):
        ax.plot(T,Y[:,j + 4], color=Col_l[j,:],linewidth=3)
        ax.plot(T,S[:,j + 4],'--', color=Col_d[j,:],linewidth=3)

    
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1e-5, Tmax)
    ax.set_ylim(1e-20, 1)
    ax.legend(label=['S5','S6','S7','P'], loc='best')
    fig.show()


elif pyruvate == 'infinite':
        
    fig, ax = plt.subplots()
    for j in range(4):
        ax.plot(T,Y[:,j], color=Col_l[j,:],linewidth=3)
        ax.plot(T,S[:,j],'--', color=Col_d[j,:],linewidth=3)
    
    ax.set_xlim(1e-2, Tmax)
    ax.set_ylim(1e-10, 1e2)
    ax.legend(label=['S1','S2','S3','S4'], loc='best')
    fig.show()

    fig, ax = plt.subplots()
    for j in range(4):
        ax.plot(T,Y[:,j + 4], color=Col_l[j,:],linewidth=3)
        ax.plot(T,S[:,j + 4],'--', color=Col_d[j,:],linewidth=3)

    ax.set_xlim(1e-2, Tmax)
    ax.set_ylim(1e-10, 1e2)
    plt.legend(label=['S5','S6','S7','P'], loc='best')
    fig.show()
else:
    raise TypeError


# switch pyruvate
#     case 'finite'
#         J_tot = 2
#     case 'infinite'
#         J_tot = 4
# end

type = 'dim'

for j in range(2):
    fig, ax = plt.subplots()
    if type == 'dim':
            ax.set_xlabel('Time (s)') #set(obj,'Interpreter','LaTex')
            ax.set_ylabel('Concentration (M)') #set(obj,'Interpreter','LaTex','Rotation',90)
    elif type == 'non-dim':
        ax.set_xlabel('t')# set(obj,'Interpreter','LaTex')
        ax.set_ylabel('Concentration') #set(obj,'Interpreter','LaTex','Rotation',90)
    else:
        raise TypeError
  

# Original matlab plotting settings:
# set(gca,'FontSize',18)
# set(findall(gcf,'type','text'),'FontSize',18)
# set(gcf, color=[1,1,1])

plt.rcParams.update({'axes.titlesize': 'large'})