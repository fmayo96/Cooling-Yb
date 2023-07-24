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
    
    #Transition rates:
    dim = 7
    b = np.zeros(dim**2, dtype=complex)
    M = np.zeros([dim**2, dim**2], dtype = complex)
    ks = [k_12, k_23, k_34, Omega_4/np.sqrt(N(beta, w)), k_56, k_67]
    ws = [w1, w2, w3, w, w4, w5]
    for i in range(1, dim-1):
        for j in range(1, dim-1):
            M[dim*i + j, dim*i+j+1] = 1j*np.sqrt(N(beta, ws[j]))*ks[j]
            M[dim*i + j, dim*i+j-1] = 1j*np.sqrt(N(beta, ws[j-1]))*ks[j-1]
            M[dim*i + j, dim*(i+1)+j] = -1j*np.sqrt(N(beta, ws[i]))*ks[i]
            M[dim*i + j, dim*(i-1)+j] = -1j*np.sqrt(N(beta, ws[i-1]))*ks[i-1]
    #===i=0
    for j in range(1, dim-1):
        M[j, dim*0+j+1] = 1j*np.sqrt(N(beta, ws[j]))*ks[j]
        M[j, dim*0+j-1] = 1j*np.sqrt(N(beta, ws[j-1]))*ks[j-1]
        M[j, dim*(0+1)+j] = -1j*np.sqrt(N(beta, ws[0]))*ks[0]
    #===i=0, j=0
    M[0,1] = 1j*np.sqrt(N(beta, ws[0]))*ks[0]
    M[0,dim] = -1j*np.sqrt(N(beta, ws[1]))*ks[1]
    #===i=dim-1
    for j in range(1, dim-1):
        M[dim*(dim-1) + j, dim*(dim-1)+j+1] = 1j*np.sqrt(N(beta, ws[j]))*ks[j]
        M[dim*(dim-1) + j, dim*(dim-1)+j-1] = 1j*np.sqrt(N(beta, ws[j-1]))*ks[j-1]
        M[dim*(dim-1) + j, dim*((dim-1)-1)+j] = -1j*np.sqrt(N(beta, ws[(dim-1)-1]))*ks[(dim-1)-1]
    #===i=dim-1, j=dim-1
    M[dim*(dim-1) + (dim-1), dim*(dim-1)+(dim-1)-1] = 1j*np.sqrt(N(beta, ws[(dim-1)-1]))*ks[(dim-1)-1]
    M[dim*(dim-1) + (dim-1), dim*((dim-1)-1)+(dim-1)] = -1j*np.sqrt(N(beta, ws[(dim-1)-1]))*ks[(dim-1)-1]   
    #===j=0
    for i in range(1, dim-1):
        M[dim*i, dim*i+1] = 1j*np.sqrt(N(beta, ws[0]))*ks[0]
        M[dim*i, dim*(i+1)] = -1j*np.sqrt(N(beta, ws[i]))*ks[i]
        M[dim*i, dim*(i-1)] = -1j*np.sqrt(N(beta, ws[i-1]))*ks[i-1]
    #===j=dim-1
    for i in range(1, dim-1):
        M[dim*i + (dim-1), dim*i+(dim-1)-1] = 1j*np.sqrt(N(beta, ws[(dim-1)-1]))*ks[(dim-1)-1]
        M[dim*i + (dim-1), dim*(i+1)+(dim-1)] = -1j*np.sqrt(N(beta, ws[i]))*ks[i]
        M[dim*i + (dim-1), dim*(i-1)+(dim-1)] = -1j*np.sqrt(N(beta, ws[i-1]))*ks[i-1]
    #===i=0, j=dim-1
    M[dim-1, dim-2] = 1j*np.sqrt(N(beta, ws[dim-2]))*ks[dim-2]
    M[dim-1, dim+dim-1] = -1j*np.sqrt(N(beta, ws[0]))*ks[0]
    #===i=dim-1, j=0
    M[dim*(dim-1), dim*(dim-1)+1] = 1j*np.sqrt(N(beta, ws[0]))*ks[0]
    M[dim*(dim-1), dim*(dim-2)] = -1j*np.sqrt(N(beta, ws[dim-2]))*ks[dim-2]
    
    for i in range(0,3):
        for j in range(4,dim):
            M[dim*i+j, dim*i+j] -= 0.5*gamma_nr

    for i in range(0, dim-1):
        M[dim*i + i + 1, dim*i + i + 1] -= 0.5*gamma_nr
    
    for i in range(0, dim-1):
        if i !=3:
            M[dim*i+i, dim*(i+1) + i+1] += gamma_nr 
    for i in range(1, dim):
        if i != 4:
            M[dim*i+i, dim*i + i] -= gamma_nr

    for i in range(4,dim):
        M[dim*i + i, dim*i + i] -= 4*gamma
    for i in range(0,4):
        for j in range(4, dim):
            M[dim*i+i, dim*j+j] += gamma

   

    print(np.linalg.det(M))


    
    return 0 

Steady_state(300,0,0.6)


