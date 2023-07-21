#Solve the equations to find the steady state 

import numpy as np 
import scipy.linalg
from parameters import *
from aux_func import *


def Steady_state(T, j_0_3, j_0_4):
    """
    This function calculates the steady state of the system:
        inputs: -T: Temperature in Kelvin.
                -j_0_3: pump intensity of the laser that drives the 3-5 transition, in MW/cm^2.
                -j_0_4: pump intensity of the laser that drives the 4-5 transition, in MW/cm^2.
        output: steady state of the system. Format =  array.
    """
    beta = 1/(kB*T)
    Omega_3 = Rabifreq(j_0_3*1000000)
    Omega_4 = Rabifreq(j_0_4*1000000)
    print('Omega_4=',Omega_4*1e-9,'ns^-1')
    #Transition rates:
    gamma_mn_1 = gamma_nr*(N(beta,w1) + 1)
    gamma_mn_2 = gamma_nr*(N(beta,w2) + 1)
    gamma_mn_3 = gamma_nr*(N(beta,w3) + 1)
    gamma_mn_5 = gamma_nr*(N(beta,w4) + 1)
    gamma_mn_6 = gamma_nr*(N(beta,w5) + 1)
    gamma_pl_1 = gamma_nr*N(beta,w1)
    gamma_pl_2 = gamma_nr*N(beta,w2)
    gamma_pl_3 = gamma_nr*N(beta,w3)
    gamma_pl_5 = gamma_nr*N(beta,w4)
    gamma_pl_6 = gamma_nr*N(beta,w5)
    b = np.zeros(12, dtype=complex)
    b[0:5] = -gamma
    b[5] = -gamma_mn_6

    M = np.zeros([12,12], dtype = complex)
    M[0,0],M[0,1],M[0,2],M[0,3] = -(gamma_pl_1 + gamma), gamma_mn_1 - gamma, -gamma, -gamma
    M[1,0],M[1,1],M[1,2],M[1,3] = gamma_pl_1 - gamma, -(gamma_mn_1 + gamma_pl_2 + gamma), gamma_mn_2 - gamma, -gamma
    M[2,0],M[2,1],M[2,2],M[2,3],M[2,8],M[2,9] = -gamma, gamma_pl_2 - gamma, -(gamma_mn_2 + gamma_pl_3 + gamma), gamma_mn_3 - gamma, 1j*Omega_3,-1j*Omega_3
    M[3,0],M[3,1],M[3,2],M[3,3],M[3,10],M[3,11] = -gamma,-gamma, gamma_pl_3 - gamma, -(gamma_mn_3 + gamma), 1j*Omega_4,-1j*Omega_4
    M[4,4],M[4,5],M[4,8],M[4,9],M[4,10],M[4,11] = -4*gamma-gamma_pl_5,gamma_mn_5,-1j*Omega_3,1j*Omega_3,-1j*Omega_4,1j*Omega_4
    M[5,0],M[5,1],M[5,2],M[5,3],M[5,4],M[5,5] = -gamma_mn_6,-gamma_mn_6,-gamma_mn_6,-gamma_mn_6,gamma_pl_5 - gamma_mn_6, -(gamma_mn_5 + gamma_mn_6 + gamma_pl_6 +4*gamma)
    M[6,6],M[6,8],M[6,11] = -1j*(E2-E3)-0.5*(gamma_mn_2+gamma_mn_3+gamma_pl_3),1j*Omega_4,-1j*Omega_3
    M[7,7],M[7,9],M[7,10] = -1j*(E3-E2)-0.5*(gamma_mn_2+gamma_mn_3+gamma_pl_3),-1j*Omega_4,1j*Omega_3
    M[8,2],M[8,4],M[8,6],M[8,8] = 1j*Omega_3,-1j*Omega_3,1j*Omega_4,-1j*(E2-E4)-0.5*(gamma_mn_2+gamma_pl_3+gamma_pl_5)
    M[9,2],M[9,4],M[9,7],M[9,9] = -1j*Omega_3,1j*Omega_3,-1j*Omega_4,1j*(E2-E4)-0.5*(gamma_mn_2+gamma_pl_3+gamma_pl_5)
    M[10,3],M[10,4],M[10,7],M[10,10] = 1j*Omega_4,-1j*Omega_4,1j*Omega_3,-1j*(E3-E4)-0.5*(gamma_mn_3-gamma_pl_5)
    M[11,3],M[11,4],M[11,6],M[11,11] = -1j*Omega_4,1j*Omega_4,-1j*Omega_3,1j*(E3-E4)-0.5*(gamma_mn_3-gamma_pl_5)


    rho_ss = scipy.linalg.solve(M,b)


    return np.append( np.real(rho_ss[0:6]), 1 - np.sum(np.real(rho_ss[0:6])) )



