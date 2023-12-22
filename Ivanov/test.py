import numpy as np 
import numpy as np
from parameters import *
from aux_func import *
import matplotlib.pyplot as plt 

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

rho = np.loadtxt('Test_j=0/ss_j=00_300.txt')

print(rho)