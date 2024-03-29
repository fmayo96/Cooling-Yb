import numpy as np 
import numpy as np
from parameters import *
from aux_func import *
import matplotlib.pyplot as plt 

TT ='150'
rho_00 = np.loadtxt(TT+'K/ss_j0=00_'+TT+'.txt')
rho_01 = np.loadtxt(TT+'K/ss_j0=01_'+TT+'.txt')
rho_02 = np.loadtxt(TT+'K/ss_j0=02_'+TT+'.txt')
rho_03 = np.loadtxt(TT+'K/ss_j0=03_'+TT+'.txt')
rho_04 = np.loadtxt(TT+'K/ss_j0=04_'+TT+'.txt')
rho_05 = np.loadtxt(TT+'K/ss_j0=05_'+TT+'.txt')
rho_06 = np.loadtxt(TT+'K/ss_j0=06_'+TT+'.txt')
rho_07 = np.loadtxt(TT+'K/ss_j0=07_'+TT+'.txt')
rho_08 = np.loadtxt(TT+'K/ss_j0=08_'+TT+'.txt')
rho_09 = np.loadtxt(TT+'K/ss_j0=09_'+TT+'.txt')
rho_10 = np.loadtxt(TT+'K/ss_j0=10_'+TT+'.txt')
rho_11 = np.loadtxt(TT+'K/ss_j0=11_'+TT+'.txt')
rho_12 = np.loadtxt(TT+'K/ss_j0=12_'+TT+'.txt')
rho_13 = np.loadtxt(TT+'K/ss_j0=13_'+TT+'.txt')
rho_14 = np.loadtxt(TT+'K/ss_j0=14_'+TT+'.txt')
rho_15 = np.loadtxt(TT+'K/ss_j0=15_'+TT+'.txt')

rho_secondlaser_00 = np.loadtxt(TT+'K_secondlaser/ss_j0=00_'+TT+'.txt')
rho_secondlaser_01 = np.loadtxt(TT+'K_secondlaser/ss_j0=01_'+TT+'.txt')
rho_secondlaser_02 = np.loadtxt(TT+'K_secondlaser/ss_j0=02_'+TT+'.txt')
rho_secondlaser_03 = np.loadtxt(TT+'K_secondlaser/ss_j0=03_'+TT+'.txt')
rho_secondlaser_04 = np.loadtxt(TT+'K_secondlaser/ss_j0=04_'+TT+'.txt')
rho_secondlaser_05 = np.loadtxt(TT+'K_secondlaser/ss_j0=05_'+TT+'.txt')
rho_secondlaser_06 = np.loadtxt(TT+'K_secondlaser/ss_j0=06_'+TT+'.txt')
rho_secondlaser_07 = np.loadtxt(TT+'K_secondlaser/ss_j0=07_'+TT+'.txt')
rho_secondlaser_08 = np.loadtxt(TT+'K_secondlaser/ss_j0=08_'+TT+'.txt')
rho_secondlaser_09 = np.loadtxt(TT+'K_secondlaser/ss_j0=09_'+TT+'.txt')
rho_secondlaser_10 = np.loadtxt(TT+'K_secondlaser/ss_j0=10_'+TT+'.txt')
rho_secondlaser_11 = np.loadtxt(TT+'K_secondlaser/ss_j0=11_'+TT+'.txt')
rho_secondlaser_12 = np.loadtxt(TT+'K_secondlaser/ss_j0=12_'+TT+'.txt')
rho_secondlaser_13 = np.loadtxt(TT+'K_secondlaser/ss_j0=13_'+TT+'.txt')
rho_secondlaser_14 = np.loadtxt(TT+'K_secondlaser/ss_j0=14_'+TT+'.txt')
rho_secondlaser_15 = np.loadtxt(TT+'K_secondlaser/ss_j0=15_'+TT+'.txt')


rho = [rho_00, rho_01, rho_02, rho_03, rho_04, rho_05, rho_06, rho_07, rho_08, rho_09, rho_10, 
       rho_11, rho_12, rho_13, rho_14, rho_15]

rho_secondlaser = [rho_secondlaser_00, rho_secondlaser_01, rho_secondlaser_02, rho_secondlaser_03, rho_secondlaser_04, rho_secondlaser_05, rho_secondlaser_06, rho_secondlaser_07, rho_secondlaser_08, rho_secondlaser_09, rho_secondlaser_10, 
       rho_secondlaser_11, rho_secondlaser_12, rho_secondlaser_13, rho_secondlaser_14, rho_secondlaser_15]



def plot_pow(rho, T, j_0_3, j_0_4, laser = 'first'):
    if laser == 'first':
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

        eps_vib_C1 = rho_0 * E3 + np.sum( E_gs * Delta_rho_gs) # J

        eps_vib_C2 = np.sum( E_es * steady_state[4:] ) + rho_0 * (E3 - E4 - .25 * np.sum(E_gs))

        n_ion = 1.3e20 # cm-3
        N_e = 13
        eta_e = 0.92 # extraction efficiency

        # Cooling power
        P_cool = eta_e * n_ion * N_e * gamma * (eps_vib_C1 + eps_vib_C2) 
        


        eps_vib_H1 = rho_0 * np.sum( .25 * E_gs ) + np.sum( E_gs * Delta_rho_gs  ) 
        eps_vib_H2 = rho_0 * (E4 - E3 + eps_vib_H1 + eps_vib_C1 )
        W_NR = 1.45  # s-1

        P_nr = n_ion * N_e *  ( gamma * eps_vib_H1 + W_NR * eps_vib_H2 ) 
        

        P_im = (j_0_3 + j_0_4) * alpha_imp *1e6

        
        P_abs = (j_0_3 + j_0_4) * (alpha_imp + alpha_rad) * 1e6
        
        P_heat = P_nr + P_im
        P_net = P_heat - P_cool
        eta = -P_net/P_abs
    
    else:
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

        eps_vib_C1 = rho_0 * E2 + np.sum( E_gs * Delta_rho_gs) # J

        eps_vib_C2 = np.sum( E_es * steady_state[4:] ) + rho_0 * (E2 - E4 - .25 * np.sum(E_gs))

        n_ion = 1.3e20 # cm-3
        N_e = 13
        eta_e = 0.92 # extraction efficiency

        # Cooling power
        P_cool = eta_e * n_ion * N_e * gamma * (eps_vib_C1 + eps_vib_C2) 
        


        eps_vib_H1 = rho_0 * np.sum( .25 * E_gs ) + np.sum( E_gs * Delta_rho_gs  ) 
        eps_vib_H2 = rho_0 * (E4 - E2 + eps_vib_H1 + eps_vib_C1 )
        W_NR = 1.45  # s-1

        P_nr = n_ion * N_e *  ( gamma * eps_vib_H1 + W_NR * eps_vib_H2 ) 
        

        P_im = (j_0_3 + j_0_4) * alpha_imp *1e6

        
        P_abs = (j_0_3 + j_0_4) * (alpha_imp + alpha_rad) * 1e6
        
        P_heat = P_nr + P_im
        P_net = P_heat - P_cool
        eta = -P_net/P_abs


    return P_net

pows = []
pows_secondlaser = []
js = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]

for i in range(0,len(js)):
    pows.append(plot_pow(rho[i], 150, 0, js[i], laser='first'))
    pows_secondlaser.append(plot_pow(rho_secondlaser[i], 150, js[i], 0, laser='second'))
    
pows = np.array(pows)
pows_secondlaser = np.array(pows_secondlaser)
print(pows)
plt.figure()
plt.plot(js,-pows, '.-C0', linewidth = 2, label = '150 K')
plt.plot(js,-pows_secondlaser, '.-C3', linewidth = 2, label = '150 K')
plt.legend()
plt.grid()
#plt.ylim([-20,400])
plt.xlabel('Pump intensity (MW/cm**2)')
plt.ylabel('Net Cooling Power')
plt.show()


Es = np.array([E0,E1,E2,E3,E4,E5,E6])

plt.figure()
plt.plot(Es, rho_00, '.')
plt.yscale('log')
plt.show()