# Calculates all components of heating and cooling power

import numpy as np
from parameters import *
from aux_func import *
from steady_state import *


def Net_Power(T, j_0_3, j_0_4):
    beta = 1/(kB*T)
    #steady_state = Steady_state(T,j_0_3,j_0_4)
    steady_state_0 = np.array([0.57955245, 0.18618485, 0.095937,   0.05952767, 0.04397912, 0.02347037,
 0.01134855])
    steady_state = np.array([0.56632714, 0.18193613, 0.09374773, 0.05816925, 0.05571185, 0.02973179,
 0.01437612])
    print(steady_state.sum())
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
    eta_e = 0.92 # extraction efficiency

    # Cooling power
    P_cool = eta_e * n_ion * N_e * rho_0 * gamma * (eps_vib_C1 + eps_vib_C2) 
    


    eps_vib_H1 = np.sum( .25 * E_gs ) + np.sum( E_gs * Delta_rho_gs / rho_0 ) 
    eps_vib_H2 = E4 - E3 + eps_vib_H1 + eps_vib_C1 
    W_NR = 1.45  # s-1

    P_nr = n_ion * N_e * rho_0 * ( gamma * eps_vib_H1 + W_NR * eps_vib_H2 ) 
    

    P_im = (j_0_3 + j_0_4) * alpha_imp *1000000

    
    P_abs = (j_0_3 + j_0_4) * (alpha_imp + alpha_rad) * 1000000
    
    P_heat = P_nr + P_im
    P_net = P_heat - P_cool
    eta = -P_net/P_abs
    
    print('Pump intensity 4 =', j_0_4, ' MW / cm2 \n')
    print('Pump intensity 3 =', j_0_3, ' MW / cm2 \n')
    print('P_nr =', P_nr, 'W / cm3')
    print('P_im =', P_im, 'W / cm3 \n')

    print('P_heat =', P_heat, 'W / cm3')
    print('P_cool =', P_cool, 'W / cm3 \n')

    print('P_net = ', P_net , 'W / cm3 \n')

    print('P_abs =', P_abs, 'W / cm3 \n')
    print('eta_c = - P_net / P_abs =', eta)
    
    return P_net

Pow_300_06 = Net_Power(300, 0, 0.6)



print(Pow_300_06)
