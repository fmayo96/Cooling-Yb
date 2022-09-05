from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 

#---Defino distribucion de Planck para los rates de disipacion---
def N_planck(beta, omega):
    return 1/(np.exp(beta*omega) - 1)

#--------------------------Parametros-----------------------------
N_steps = 100000
dt = 1e-3
tf = N_steps*dt
w = 1
w_laser = 10*w 
g = 0.5
gamma = 0.1
alpha = np.linspace(0,10,10)
beta = 10

gamma_pl = g**2*N_planck(beta,w)
gamma_mn = g**2*(N_planck(beta,w) + 1)

heat_current = []
heat_current_2 = []
#----------One Laser--------------------
for j in range(len(alpha)):

    rho = np.zeros([N_steps,7])
    rho_dot = np.zeros(7)
    rho_thermal = np.zeros(7)
    for i in range(0,4):
        rho_thermal[i] = np.exp(-beta*i*w)
    rho_thermal[4] = np.exp(-beta*(3*w+w_laser))
    rho_thermal[5] = np.exp(-beta*(4*w+w_laser))
    rho_thermal[6] = np.exp(-beta*(5*w+w_laser))
    rho_thermal /= np.sum(rho_thermal)
    rho[0] = rho_thermal
    M = np.array([[-gamma_pl,gamma_mn,0,0,gamma,gamma,gamma],[gamma_pl,-(gamma_pl+gamma_mn),gamma_mn,0,gamma,gamma,gamma],[0,gamma_pl,-(gamma_pl + gamma_mn),gamma_mn,gamma,gamma,gamma],[0,0,gamma_pl,-(alpha[j] + gamma_mn),gamma,gamma,gamma],[0,0,0,alpha[j],-(4*gamma + gamma_pl),gamma_mn,0],[0,0,0,0,gamma_pl,-(4*gamma + gamma_pl + gamma_mn),gamma_mn],[0,0,0,0,0,gamma_pl,-(4*gamma + gamma_mn)]])

    for i in range(N_steps-1):
        rho[i+1] = rho[i] + dt*np.dot(M,rho[i])
    t = np.linspace(0,tf,N_steps)
    x = [0,w,2*w,3*w,3*w+w_laser,4*w+w_laser,5*w+w_laser]
    #------------------------------Heat-----------------------------
    E1 = 0
    E2 = w
    E3 = 2*w 
    E4 = 3*w 
    E5 = 3*w + w_laser
    E6 = 4*w + w_laser
    E7 = 5*w + w_laser

    M_heat = np.array([[-gamma_pl*E1,gamma_mn*E1,0,0,0,0,0],[gamma_pl*E2,-(gamma_pl*E2+gamma_mn*E2),gamma_mn*E2,0,0,0,0],[0,gamma_pl*E3,-E3*(gamma_pl + gamma_mn),E3*gamma_mn,0,0,0],[0,0,gamma_pl*E4,-(gamma_mn)*E4,0,0,0],[0,0,0,0,-(gamma_pl)*E5,gamma_mn*E5,0],[0,0,0,0,gamma_pl*E6,-(gamma_pl + gamma_mn)*E6,gamma_mn*E6],[0,0,0,0,0,gamma_pl*E7,-(gamma_mn)*E7]])

    Q_i = np.sum(np.dot(M_heat,rho[0]))
    Q_f = np.sum(np.dot(M_heat,rho[-1]))
    heat_current.append(Q_f)

#----------Two Laser--------------------
for j in range(len(alpha)):

    rho = np.zeros([N_steps,7])
    rho_dot = np.zeros(7)
    rho_thermal = np.zeros(7)
    for i in range(0,4):
        rho_thermal[i] = np.exp(-beta*i*w)
    rho_thermal[4] = np.exp(-beta*(3*w+w_laser))
    rho_thermal[5] = np.exp(-beta*(4*w+w_laser))
    rho_thermal[6] = np.exp(-beta*(5*w+w_laser))
    rho_thermal /= np.sum(rho_thermal)
    rho[0] = rho_thermal
    M = np.array([[-gamma_pl,gamma_mn,0,0,gamma,gamma,gamma],[gamma_pl,-(gamma_pl+gamma_mn),gamma_mn,0,gamma,gamma,gamma],[0,gamma_pl,-(gamma_pl + gamma_mn + alpha[j]),gamma_mn,gamma,gamma,gamma],[0,0,gamma_pl,-(alpha[j] + gamma_mn),gamma,gamma,gamma],[0,0,alpha[j],alpha[j],-(4*gamma + gamma_pl),gamma_mn,0],[0,0,0,0,gamma_pl,-(4*gamma + gamma_pl + gamma_mn),gamma_mn],[0,0,0,0,0,gamma_pl,-(4*gamma + gamma_mn)]])

    for i in range(N_steps-1):
        rho[i+1] = rho[i] + dt*np.dot(M,rho[i])
    t = np.linspace(0,tf,N_steps)
    x = [0,w,2*w,3*w,3*w+w_laser,4*w+w_laser,5*w+w_laser]
    #------------------------------Heat-----------------------------
    E1 = 0
    E2 = w
    E3 = 2*w 
    E4 = 3*w 
    E5 = 3*w + w_laser
    E6 = 4*w + w_laser
    E7 = 5*w + w_laser

    M_heat = np.array([[-gamma_pl*E1,gamma_mn*E1,0,0,0,0,0],[gamma_pl*E2,-(gamma_pl*E2+gamma_mn*E2),gamma_mn*E2,0,0,0,0],[0,gamma_pl*E3,-E3*(gamma_pl + gamma_mn),E3*gamma_mn,0,0,0],[0,0,gamma_pl*E4,-(gamma_mn)*E4,0,0,0],[0,0,0,0,-(gamma_pl)*E5,gamma_mn*E5,0],[0,0,0,0,gamma_pl*E6,-(gamma_pl + gamma_mn)*E6,gamma_mn*E6],[0,0,0,0,0,gamma_pl*E7,-(gamma_mn)*E7]])

    Q_i = np.sum(np.dot(M_heat,rho[0]))
    Q_f = np.sum(np.dot(M_heat,rho[-1]))
    heat_current_2.append(Q_f)



plt.figure()
plt.plot(alpha,heat_current,'o', label = '1 Laser')
plt.plot(alpha,heat_current_2,'C3o', label = '2 Laser')
plt.xlabel(r'$\alpha$', fontsize = 12)
plt.ylabel(r'$\dot{Q}$', fontsize = 12)
plt.legend(fontsize = 11)
plt.show()



"""
    plt.figure()
    plt.plot(x,rho[0],linewidth = 2, label = r'$\rho_i$')
    plt.plot(x,rho[-1],'o',linewidth = 2, label = r'$\rho_f$')
    plt.xlabel('Energy', fontsize = 14)
    plt.ylabel('Populations', fontsize = 14)
    plt.yscale('log')
    plt.legend(fontsize = 12)
    plt.show()
    
    plt.figure()
    plt.plot(t,rho[:,0], linewidth = 2, label = r'$\rho_1$')
    plt.plot(t,rho[:,1], linewidth = 2, label = r'$\rho_2$')
    plt.plot(t,rho[:,2], linewidth = 2, label = r'$\rho_3$')
    plt.plot(t,rho[:,3], linewidth = 2, label = r'$\rho_4$')
    plt.plot(t,rho[:,4], linewidth = 2, label = r'$\rho_5$')
    plt.plot(t,rho[:,5], linewidth = 2, label = r'$\rho_6$')
    plt.plot(t,rho[:,6], linewidth = 2, label = r'$\rho_7$')
    plt.xlabel('Time', fontsize = 14)
    plt.legend(fontsize = 12)
    plt.show()
"""