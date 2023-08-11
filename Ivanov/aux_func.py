#All short functions are defined in this file

import numpy as np 
from parameters import *
from numba import njit 

@njit
def Rabifreq(j_0):
    E_0 = np.sqrt(2*j_0/(c*n*eps0))
    return d*E_0/hbar
@njit
def N(beta,omega):
    return 1/(np.exp(beta*omega)-1)
@njit
def Thermal_state(beta, H):
    return np.exp(-beta*H)/np.sum(np.exp(-beta*H))

