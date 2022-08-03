import numpy as np
from datetime import time
import matplotlib.pyplot as plt

[T,Y,k,ep,pyruvate,Tmax] = Full_dynamics

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
[S,t_long,S_long,S1_med] = asymptotic_solutions(T,k,ep)
print(time.time-start)

# size(S)

# S1 = S(1,:)
# S2 = S(2,:)
# S3 = S(3,:)
# S4 = S(4,:)
​
if pyruvate == 'finite':
        fig, ax = plt.subplots(1)
        for j in range(1,4):
            ax.plot(T,Y(:,j),'Color',Col_l(j,:),'LineWidth',3, label=f'S{j}')
            ax.set_yscale('log')
            ax.set_xscale('log')

        ax.set_xlim(1e-5, Tmax)
        ax.set_ylim(1e-20, 1)
        plt.legend(loc='best')
        
        for j in range(1,4):
            ax.plot(T,S[:,j],'--','Color',Col_d[j,:],'LineWidth',3)
        
        
        fig, ax = plt.subplots(2)
        for j = 1:4
        loglog(T,Y(:,j + 4),'Color',Col_l(j,:),'LineWidth',3)
        hold on
        end
        
        for j = 1:4
        plt.plot(T,S(:,j + 4),'--','Color',Col_d(j,:),'LineWidth',3)
        end
        
        axis([1e-5 Tmax 1e-20 1])
        plt.legend('S5','S6','S7','P')
        
        
    case 'infinite'
        
        
        fig, ax = plt.subplots(1)
        for j = 1:4
        loglog(T,Y(:,j),'Color',Col_l(j,:),'LineWidth',3)
        hold on
        end
        
        for j = 1:4
        plt.plot(T,S(:,j),'--','Color',Col_d(j,:),'LineWidth',3)
        end
        
        
        axis([1e-2 Tmax 1e-10 1e2])
        plt.legend('S1','S2','S3','S4')
        
        fig, ax = plt.subplots(2)
        for j = 1:4
        loglog(T,Y(:,j + 4),'Color',Col_l(j,:),'LineWidth',3)
        hold on
        end
        
        for j = 1:4
        plt.plot(T,S(:,j + 4),'--','Color',Col_d(j,:),'LineWidth',3)
        end
        axis([1e-2 Tmax 1e-10 1e2])
        plt.legend('S5','S6','S7','P')
        
        
    otherwise
        error('Unexpected type. Cannot proceed')
end
​
​
​
# switch pyruvate
#     case 'finite'
#         J_tot = 2
#     case 'infinite'
#         J_tot = 4
# end

type = 'dim'

for j = 1:2
    fig, ax = plt.subplots(j)
switch type
    case 'dim'
        obj = xlabel('Time (s)') #set(obj,'Interpreter','LaTex')
        obj = ylabel('Concentration (M)') #set(obj,'Interpreter','LaTex','Rotation',90)
    case 'non-dim'
        obj = xlabel('t')# set(obj,'Interpreter','LaTex')
        obj = ylabel('Concentration') #set(obj,'Interpreter','LaTex','Rotation',90)
    otherwise
        error('Unexpected type. Cannot proceed')
end

set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
set(gcf,'Color',[1,1,1])
end