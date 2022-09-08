import numpy as np 
import matplotlib.pyplot as plt 


#---Defino distribucion de Planck para los rates de disipacion---
def N_planck(beta, omega):
    return 1/(np.exp(beta*omega) - 1)

#--------------------------Parametros-----------------------------
N_steps = 10000
dt = 1e-3
tf = N_steps*dt
w = 1
w_laser = 10*w 
g = 0.5
gamma = 0.1
alpha = 0.5
#T = np.linspace(5e-2, 1, 20)
beta = np.linspace(1,20,20)
T = 1/beta


heat_current = []
heat_current_2 = []
#----------One Laser--------------------
for j in range(len(beta)):
    gamma_pl = g**2*N_planck(beta[j],w)
    gamma_mn = g**2*(N_planck(beta[j],w) + 1)
    rho = np.zeros([N_steps,7])
    rho_dot = np.zeros(7)
    rho_thermal = np.zeros(7)
    for i in range(0,4):
        rho_thermal[i] = np.exp(-beta[j]*i*w)
    rho_thermal[4] = np.exp(-beta[j]*(3*w+w_laser))
    rho_thermal[5] = np.exp(-beta[j]*(4*w+w_laser))
    rho_thermal[6] = np.exp(-beta[j]*(5*w+w_laser))
    rho_thermal /= np.sum(rho_thermal)
    rho[0] = rho_thermal
    M = np.array([[-gamma_pl,gamma_mn,0,0,gamma,gamma,gamma],[gamma_pl,-(gamma_pl+gamma_mn),gamma_mn,0,gamma,gamma,gamma],[0,gamma_pl,-(gamma_pl + gamma_mn),gamma_mn,gamma,gamma,gamma],[0,0,gamma_pl,-(alpha + gamma_mn),gamma,gamma,gamma],[0,0,0,alpha,-(4*gamma + gamma_pl),gamma_mn,0],[0,0,0,0,gamma_pl,-(4*gamma + gamma_pl + gamma_mn),gamma_mn],[0,0,0,0,0,gamma_pl,-(4*gamma + gamma_mn)]])

    for i in range(N_steps-1):
        K1 = np.dot(M,rho[i])
        K2 = np.dot(M,rho[i] + dt*K1/2)
        K3 = np.dot(M,rho[i] + dt*K2/2)
        K4 = np.dot(M,rho[i] + dt*K3)
        rho[i+1] = rho[i] + dt/6*(K1 + 2*K2 + 2*K3 + K4)
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
for j in range(len(beta)):
    gamma_pl = g**2*N_planck(beta[j],w)
    gamma_mn = g**2*(N_planck(beta[j],w) + 1)
    rho = np.zeros([N_steps,7])
    rho_dot = np.zeros(7)
    rho_thermal = np.zeros(7)
    for i in range(0,4):
        rho_thermal[i] = np.exp(-beta[j]*i*w)
    rho_thermal[4] = np.exp(-beta[j]*(3*w+w_laser))
    rho_thermal[5] = np.exp(-beta[j]*(4*w+w_laser))
    rho_thermal[6] = np.exp(-beta[j]*(5*w+w_laser))
    rho_thermal /= np.sum(rho_thermal)
    rho[0] = rho_thermal
    M = np.array([[-gamma_pl,gamma_mn,0,0,gamma,gamma,gamma],[gamma_pl,-(gamma_pl+gamma_mn),gamma_mn,0,gamma,gamma,gamma],[0,gamma_pl,-(gamma_pl + gamma_mn + alpha),gamma_mn,gamma,gamma,gamma],[0,0,gamma_pl,-(alpha + gamma_mn),gamma,gamma,gamma],[0,0,alpha,alpha,-(4*gamma + gamma_pl),gamma_mn,0],[0,0,0,0,gamma_pl,-(4*gamma + gamma_pl + gamma_mn),gamma_mn],[0,0,0,0,0,gamma_pl,-(4*gamma + gamma_mn)]])

    for i in range(N_steps-1):
        K1 = np.dot(M,rho[i])
        K2 = np.dot(M,rho[i] + dt*K1/2)
        K3 = np.dot(M,rho[i] + dt*K2/2)
        K4 = np.dot(M,rho[i] + dt*K3)
        rho[i+1] = rho[i] + dt/6*(K1 + 2*K2 + 2*K3 + K4)
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
plt.plot(T,heat_current,'o', label = '1 Laser')
plt.plot(T,heat_current_2,'C3o', label = '2 Laser')
plt.xlabel(r'$T_{crystal}$', fontsize = 12)
plt.ylabel(r'$\dot{Q}$', fontsize = 12)
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 11)
plt.show()

