import numpy as np 
import numpy as np
from parameters import *
from aux_func import *
import matplotlib.pyplot as plt 

TT ='300'

def read_file(file_path):
    with open(file_path, "r") as f:
        for line in f:
            rhoString = line.strip().split(",")
        output = [float(num) for num in rhoString]
    return output

rho_00_300 = read_file("outs/results_temp=300.0_j0=0.00.txt")
rho_01_300 = read_file("outs/results_temp=300.0_j0=0.10.txt")
rho_02_300 = read_file("outs/results_temp=300.0_j0=0.20.txt")
rho_03_300 = read_file("outs/results_temp=300.0_j0=0.30.txt")
rho_04_300 = read_file("outs/results_temp=300.0_j0=0.40.txt")
rho_05_300 = read_file("outs/results_temp=300.0_j0=0.50.txt")
rho_06_300 = read_file("outs/results_temp=300.0_j0=0.60.txt")
rho_07_300 = read_file("outs/results_temp=300.0_j0=0.70.txt")
rho_08_300 = read_file("outs/results_temp=300.0_j0=0.80.txt")
rho_09_300 = read_file("outs/results_temp=300.0_j0=0.90.txt")
rho_10_300 = read_file("outs/results_temp=300.0_j0=1.00.txt")
rho_11_300 = read_file("outs/results_temp=300.0_j0=1.10.txt")
rho_12_300 = read_file("outs/results_temp=300.0_j0=1.20.txt")
rho_13_300 = read_file("outs/results_temp=300.0_j0=1.30.txt")
rho_14_300 = read_file("outs/results_temp=300.0_j0=1.40.txt")
rho_15_300 = read_file("outs/results_temp=300.0_j0=1.50.txt")
rho_16_300 = read_file("outs/results_temp=300.0_j0=1.60.txt")
rho_17_300 = read_file("outs/results_temp=300.0_j0=1.70.txt")
rho_18_300 = read_file("outs/results_temp=300.0_j0=1.80.txt")

rho_00_300_new = read_file("new_results/results_temp=300.0_j0=0.00.txt")
rho_01_300_new = read_file("new_results/results_temp=300.0_j0=0.10.txt")
rho_02_300_new = read_file("new_results/results_temp=300.0_j0=0.20.txt")
rho_03_300_new = read_file("new_results/results_temp=300.0_j0=0.30.txt")
rho_04_300_new = read_file("new_results/results_temp=300.0_j0=0.40.txt")
rho_05_300_new = read_file("new_results/results_temp=300.0_j0=0.50.txt")
rho_06_300_new = read_file("new_results/results_temp=300.0_j0=0.60.txt")
rho_07_300_new = read_file("new_results/results_temp=300.0_j0=0.70.txt")
rho_08_300_new = read_file("new_results/results_temp=300.0_j0=0.80.txt")
rho_09_300_new = read_file("new_results/results_temp=300.0_j0=0.90.txt")
rho_10_300_new = read_file("new_results/results_temp=300.0_j0=1.00.txt")
rho_11_300_new = read_file("new_results/results_temp=300.0_j0=1.10.txt")
rho_12_300_new = read_file("new_results/results_temp=300.0_j0=1.20.txt")
rho_13_300_new = read_file("new_results/results_temp=300.0_j0=1.30.txt")
rho_14_300_new = read_file("new_results/results_temp=300.0_j0=1.40.txt")
rho_15_300_new = read_file("new_results/results_temp=300.0_j0=1.50.txt")
rho_16_300_new = read_file("new_results/results_temp=300.0_j0=1.60.txt")
rho_17_300_new = read_file("new_results/results_temp=300.0_j0=1.70.txt")
rho_18_300_new = read_file("new_results/results_temp=300.0_j0=1.80.txt")




rho_300 = [rho_00_300, rho_01_300, rho_02_300, rho_03_300, rho_04_300, rho_05_300, rho_06_300, rho_07_300, rho_08_300, rho_09_300, rho_10_300, rho_11_300, rho_12_300, rho_13_300, rho_14_300, rho_15_300, rho_16_300, rho_17_300, rho_18_300]
rho_300_new = [rho_00_300_new, rho_01_300_new, rho_02_300_new, rho_03_300_new, rho_04_300_new, rho_05_300_new, rho_06_300_new, rho_07_300_new, rho_08_300_new, rho_09_300_new, rho_10_300_new, rho_11_300_new, rho_12_300_new, rho_13_300_new, rho_14_300_new, rho_15_300_new, rho_16_300_new, rho_17_300_new, rho_18_300_new]

# rho_secondlaser = [rho_secondlaser_00, rho_secondlaser_01, rho_secondlaser_02, rho_secondlaser_03, rho_secondlaser_04, rho_secondlaser_05, rho_secondlaser_06, rho_secondlaser_07, rho_secondlaser_08, rho_secondlaser_09, rho_secondlaser_10, 
#        rho_secondlaser_11, rho_secondlaser_12, rho_secondlaser_13, rho_secondlaser_14, rho_secondlaser_15]



def plot_pow(rho, initial_state, T, j_0_3, j_0_4, laser = 'first'):
    if laser == 'first':
        beta = 1/(kB*T)
        steady_state = rho
        
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

pows_300 = []
pows_300_new = []
js = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5, 1.6, 1.7, 1.8]

for i in range(0,len(js)):
    pows_300.append(plot_pow(rho_300[i], rho_300[0], 300, 0, js[i], laser='first'))
    pows_300_new.append(plot_pow(rho_300_new[i], rho_300_new[0], 300, 0, js[i], laser='first'))
pows_300 = np.array(pows_300)
pows_300_new = np.array(pows_300_new)

plt.figure()
plt.plot(js,-pows_300, '.-C0', linewidth = 2, label = '300 K')
plt.plot(js,-pows_300_new, '.-C3', linewidth = 2, label = '300 K new')
#plt.plot(js,-pows_secondlaser, '.-C3', linewidth = 2, label = '150 K')
plt.legend()
plt.grid()
#plt.ylim([-20,400])
plt.xlabel('Pump intensity (MW/cm**2)')
plt.ylabel('Net Cooling Power')
plt.show()


Es = np.array([E0,E1,E2,E3,E4,E5,E6])

# plt.figure()
# plt.plot(Es, rho_00, '.')
# plt.yscale('log')
# plt.show()
