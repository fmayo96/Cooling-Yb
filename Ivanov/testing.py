import numpy as np 
from parameters import *
from aux_func import *

gamma_nr = 1
gamma = 2

rho = np.eye(7)
for i in range(4,7):
    rho[i,i] = 0
def nonRad(rho):
    diag = np.diag(rho)
    excitations = np.zeros(7, dtype = np.complex128)
    decay = np.zeros(7, dtype = np.complex128)
    for i in range(6):
        if i !=3:
            excitations[i] = gamma_nr*diag[i+1]
    for i in range(1,7):
        if i != 4:
            decay[i] = -gamma_nr*diag[i]
    out = np.diag(excitations) + np.diag(decay)
    return out

def spontEm(rho):
    diag = np.diag(rho)
    _spontEm = np.zeros((7,7), dtype = np.complex128)
    for i in range(4):
        _spontEm[i,i] = gamma * np.sum(diag[4:])
    for i in range(4, 7):
        _spontEm[i,i] = -4*gamma*diag[i]
    return _spontEm


beta = 1/(kB*150)
Omega_4 = 0
Omega_3 = 0
ks = [k_12, k_23, k_34, Omega_4, k_56, k_67]
ws = [w1, w2, w3, w, w4, w5]
###=========Hamiltonian=============
H_ev = np.zeros((7,7), dtype = np.complex128)
for i in range(6):
    H_ev[i, i+1] = ks[i]*(N(beta, ws[i])**0.5)
for i in range(1,7):
    H_ev[i, i-1] = ks[i-1]*(N(beta, ws[i-1])**0.5)
H_ev[3,4] = Omega_4
H_ev[4,3] = Omega_4
#H_ev[2,4] = Omega_3
#H_ev[4,2] = Omega_3

print(H_ev)

