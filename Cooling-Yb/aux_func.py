#All short functions are defined in this file

import numpy as np 
from parameters import *

def Rabifreq(j_0):
    E_0 = np.sqrt(2*j_0/(c*n*eps0))
    print('E_0=',E_0)
    return d*E_0/hbar

def N(beta,omega):
    return 1/(np.exp(beta*omega)-1)

def Thermal_state(beta, H):
    return np.diag(np.diag(np.exp(-beta*H)))/np.sum(np.diag(np.exp(-beta*H)))

