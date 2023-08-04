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
    b = np.zeros(dim**2-1, dtype = complex)
    M = np.zeros([dim**2-1, dim**2-1], dtype = complex)
    Unit = np.zeros([dim**2, dim**2], dtype= complex) #Uninatry part of the master equation
    Diss = np.zeros([dim**2, dim**2], dtype= complex) #Dissipative part of the master equation
    ks = [k_12, k_23, k_34, Omega_4/N(beta,w)**0.5, k_56, k_67]
    ws = [w1, w2, w3, w, w4, w5]
    
    
    for i in range(dim):
        for j in range(dim):
            if i>0 and j>0 and i<6 and j<6:
                Unit[dim*i+j, dim*i+j-1] = 1j*ks[j-1]*N(beta, ws[j-1])**0.5  
                Unit[dim*i+j, dim*i+j+1] = 1j*ks[j]*N(beta,ws[j])**0.5
                Unit[dim*i+j, dim*(i-1)+j] = -1j*ks[i-1]*N(beta, ws[i-1])**0.5
                Unit[dim*i+j, dim*(i+1)+j] = -1j*ks[i]*N(beta, ws[i])**0.5
            elif i>0 and j==0 and i < 6:
                Unit[dim*i+j, dim*i+j+1] = 1j*ks[j]*N(beta,ws[j])**0.5
                Unit[dim*i+j, dim*(i-1)+j] = -1j*ks[i-1]*N(beta, ws[i-1])**0.5
                Unit[dim*i+j, dim*(i+1)+j] = -1j*ks[i]*N(beta, ws[i])**0.5
            elif i==0 and j>0 and j < 6:
                Unit[dim*i+j, dim*i+j-1] = 1j*ks[j-1]*N(beta, ws[j-1])**0.5  
                Unit[dim*i+j, dim*i+j+1] = 1j*ks[j]*N(beta,ws[j])**0.5
                Unit[dim*i+j, dim*(i+1)+j] = -1j*ks[i]*N(beta, ws[i])**0.5
            elif i==0 and j==0:
                Unit[dim*i+j, dim*i+j+1] = 1j*ks[j]*N(beta,ws[j])**0.5
                Unit[dim*i+j, dim*(i+1)+j] = -1j*ks[i]*N(beta, ws[i])**0.5
            elif i == 6 and j > 0 and j < 5:
                Unit[dim*i+j, dim*i+j-1] = 1j*ks[j-1]*N(beta, ws[j-1])**0.5  
                Unit[dim*i+j, dim*i+j+1] = 1j*ks[j]*N(beta,ws[j])**0.5
                Unit[dim*i+j, dim*(i-1)+j] = -1j*ks[i-1]*N(beta, ws[i-1])**0.5
            elif i > 0 and i < 5 and j==6:
                Unit[dim*i+j, dim*i+j-1] = 1j*ks[j-1]*N(beta, ws[j-1])**0.5  
                Unit[dim*i+j, dim*(i-1)+j] = -1j*ks[i-1]*N(beta, ws[i-1])**0.5
                Unit[dim*i+j, dim*(i+1)+j] = -1j*ks[i]*N(beta, ws[i])**0.5
        Unit[dim*5+6, dim*5+5], Unit[dim*5+6, dim*4+6] = 1j*ks[5]*N(beta, ws[5])**0.5, -1j*ks[4]*N(beta, ws[4])**0.5
        Unit[dim*6+5, dim*6+4], Unit[dim*5+6, dim*5+5] = 1j*ks[4]*N(beta, ws[4])**0.5, -1j*ks[5]*N(beta, ws[5])**0.5
        for i in range(6):
            Unit[dim*5+6, dim*i+i] = 1j*ks[5]*N(beta, ws[5])**0.5
            Unit[dim*5+6, dim*i+i] = 1j*ks[5]*N(beta, ws[5])**0.5
    print(Unit[:-1, :-1])
    print(np.linalg.det(Unit[:-1, :-1]))
    return 0 

Steady_state(300,0,0.6)


