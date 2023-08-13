import numpy as np 
import numpy as np
from parameters import *
from aux_func import *
import matplotlib.pyplot as plt 

rho_01 = np.loadtxt('300K/ss_j0=01.txt')
rho_02 = np.loadtxt('300K/ss_j0=02.txt')
rho_03 = np.loadtxt('300K/ss_j0=03.txt')
rho_04 = np.loadtxt('300K/ss_j0=04.txt')
rho_05 = np.loadtxt('300K/ss_j0=05.txt')
rho_06 = np.loadtxt('300K/steadystate.txt')
rho_07 = np.loadtxt('300K/ss_j0=07.txt')
rho_08 = np.loadtxt('300K/ss_j0=08.txt')
rho_09 = np.loadtxt('300K/ss_j0=09.txt')
rho_10 = np.loadtxt('300K/ss_j0=1.txt')
rho_11 = np.loadtxt('300K/ss_j0=11.txt')
rho_12 = np.loadtxt('300K/ss_j0=12.txt')
rho_13 = np.loadtxt('300K/ss_j0=13.txt')
rho_14 = np.loadtxt('300K/ss_j0=14.txt')
rho_15 = np.loadtxt('300K/ss_j0=15.txt')


rho_01_150 = np.loadtxt('150K/ss_j0=01_150.txt')
rho_02_150 = np.loadtxt('150K/ss_j0=02_150.txt')
rho_03_150 = np.loadtxt('150K/ss_j0=03_150.txt')
rho_04_150 = np.loadtxt('150K/ss_j0=04_150.txt')
rho_05_150 = np.loadtxt('150K/ss_j0=05_150.txt')
rho_06_150 = np.loadtxt('150K/ss_j0=06_150.txt')
rho_07_150 = np.loadtxt('150K/ss_j0=07_150.txt')
rho_08_150 = np.loadtxt('150K/ss_j0=08_150.txt')
rho_09_150 = np.loadtxt('150K/ss_j0=09_150.txt')
rho_10_150 = np.loadtxt('150K/ss_j0=10_150.txt')
rho_11_150 = np.loadtxt('150K/ss_j0=11_150.txt')
rho_12_150 = np.loadtxt('150K/ss_j0=12_150.txt')
rho_13_150 = np.loadtxt('150K/ss_j0=13_150.txt')
rho_14_150 = np.loadtxt('150K/ss_j0=14_150.txt')
rho_15_150 = np.loadtxt('150K/ss_j0=15_150.txt')


rho = [rho_01, rho_02, rho_03, rho_04, rho_05, rho_06, rho_07, rho_08, rho_09, rho_10, 
       rho_11, rho_12, rho_13, rho_14, rho_15]

rho_150 = [rho_01_150, rho_02_150, rho_03_150, rho_04_150, rho_05_150, rho_06_150, rho_07_150, rho_08_150, rho_09_150, rho_10_150, 
       rho_11_150, rho_12_150, rho_13_150, rho_14_150, rho_15_150]


def plot_pow(rho, T, j_0_4):
    j_0_3 = 0
    beta = 1/(kB*T)
    steady_state = rho
    initial_state = Thermal_state(beta,H)
    
    rho_0 = np.sum(steady_state[4:])

    Delta_rho_gs = np.array([steady_state[0] - initial_state[0],
                            steady_state[1] - initial_state[1],
                            steady_state[2] - initial_state[2],
                            steady_state[3] - initial_state[3]])

    Delta_rho_es = np.array([steady_state[4] - initial_state[4],
                            steady_state[5] - initial_state[5],
                            steady_state[6] - initial_state[6]])

    eps_vib_C1 = E3 + np.sum( E_gs * Delta_rho_gs / rho_0 ) # J

    eps_vib_C2 = np.sum( E_es * steady_state[4:] / rho_0 ) + E3 - E4 - .25 * np.sum(E_gs)

    n_ion = 1.3e20 # cm-3
    N_e = 13
    eta_e = 0.95 # extraction efficiency

    # Cooling power
    P_cool = eta_e * n_ion * N_e * rho_0 * gamma * (eps_vib_C1 + eps_vib_C2) 
    


    eps_vib_H1 = np.sum( .25 * E_gs ) + np.sum( E_gs * Delta_rho_gs / rho_0 ) 
    eps_vib_H2 = E4 - E3 + eps_vib_H1 + eps_vib_C1 
    W_NR = 1.45  # s-1

    P_nr = n_ion * N_e * rho_0 * ( gamma * eps_vib_H1 + W_NR * eps_vib_H2 ) 
    

    P_im = (j_0_3 + j_0_4) * alpha_imp *1e6

    
    P_abs = (j_0_3 + j_0_4) * (alpha_imp + alpha_rad) * 1e6
    
    P_heat = P_nr + P_im
    P_net = P_heat - P_cool
    eta = -P_net/P_abs
    
    return P_net

pows = [0]
pows_150 = [0]

js = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]

for i in range(0,len(js)-1):
    pows.append(plot_pow(rho[i], 300, js[i+1]))
    pows_150.append(plot_pow(rho_150[i], 150, js[i+1]))
pows = np.array(pows)
pows_150 = np.array(pows_150)
plt.figure()
plt.plot(js,-pows, '.-C0', linewidth = 2, label = '300 K')
plt.plot(js,-pows_150, '.-C1', linewidth = 2, label = '150 K')
plt.legend()
plt.grid()
plt.ylim([0,400])
plt.xlabel('Pump intensity (MW/cm**2)')
plt.ylabel('Net Cooling Power')
plt.show()