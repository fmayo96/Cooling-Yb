#All short functions are defined in this file

import numpy as np 
from .parameters import Parameters as P


def Rabifreq(j_0, d):
    E_0 = np.sqrt(2*j_0*1e6/(P.c*P.n*P.eps0))
    return d*E_0/P.hbar

def N_therm(beta,omega):
    return 1/(np.exp(beta*omega)-1)

def Thermal_state(beta, H):
    state = np.exp(-beta*H)
    state[-3:] = [0,0,0]
    return state / np.sum(state)



def Hamiltonian(j_0_3, j_0_4):
    Omega_4 = Rabifreq(j_0_4, P.d)
    Omega_3 = Rabifreq(j_0_3, P.d2)
    H = np.zeros((7,7))
    H[3,4] = Omega_4
    H[4,3] = Omega_4
    H[2,4] = Omega_3
    H[4,2] = Omega_3
    return H

