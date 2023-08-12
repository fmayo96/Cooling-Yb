#Solve the equations to find the steady state 
import time
import numpy as np 
from parameters import *
from aux_func import *
from numba import njit

@njit
def Steady_state(T, j_0_3, j_0_4):
    """
    This function calculates the steady state of the system:
        inputs: -T: Temperature in Kelvin.
                -j_0_3: pump intensity of the laser that drives the 3-5 transition, in MW/cm^2.
                -j_0_4: pump intensity of the laser that drives the 4-5 transition, in MW/cm^2.
        output: steady state of the system. Format =  array.
    """
    beta = 1/(kB*300)
    Omega_4 = Rabifreq(j_0_4*1e6)
    ks = [k_12, k_23, k_34, Omega_4, k_56, k_67]
    ws = [w1, w2, w3, w, w4, w5]
    ###=========Hamiltonian=============
    H_ev = np.zeros((7,7), dtype = np.complex128)
    for i in range(6):
        H_ev[i, i+1] = ks[i]*N(beta, ws[i])**0.5
    for i in range(1,7):
        H_ev[i, i-1] = ks[i-1]*N(beta, ws[i-1])**0.5
    H_ev[3,4] = Omega_4
    H_ev[4,3] = Omega_4
    ###=======Decoherence Matrix========
    def decoh(rho):
        Decoh = np.zeros((7,7), dtype = np.complex128)
        for i in range(4):
            for j in range(4,7):
                Decoh[i,j] = -0.5 * gamma_nr
                Decoh[j,i] = -0.5 * gamma_nr
        for i in range(6):
            Decoh[i, i+1] = -0.5*gamma_nr
            Decoh[i+1, i] = -0.5*gamma_nr
        out = Decoh*rho
        return out
    ###======Non radiative==============
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
    #===========Spontaneous emission=======
    def spontEm(rho):
        diag = np.diag(rho)
        _spontEm = np.zeros((7,7), dtype = np.complex128)
        for i in range(4):
            _spontEm[i,i] = gamma * np.sum(diag[4:])
        for i in range(4, 7):
            _spontEm[i,i] = -4*gamma*diag[i]
        return _spontEm

    def Lliouv(rho):
        lliouv = -1j*H_ev@rho + 1j*rho@H_ev + decoh(rho) + nonRad(rho) + spontEm(rho)
        return lliouv
    tf = 1/314.15
    dt = 5e-13
    Nsteps = int(tf/dt)
    rho_i = Thermal_state(beta, H)
    print(rho_i)
    rho_i = np.diag(rho_i)
    rho = np.zeros((7,7), dtype = np.complex128)
    rho = rho_i
    for n in range(Nsteps-1):
        rho = rho + dt*Lliouv(rho)
        #if n % 1000000 == 0:
        #   print(f"n={n} | rho:{np.diag(np.real(rho))}")
    rho_ss = np.diag(np.real(rho))
    return rho_ss


start = time.time()
rho_ss = Steady_state(300,0,1.5)
print(rho_ss)
end = time.time()
print(f"Total time = {end - start}")
np.savetxt("ss_j0=15.txt", rho_ss)


