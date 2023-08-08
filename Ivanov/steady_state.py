#Solve the equations to find the steady state 

import numpy as np 
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
    Omega_3 = Rabifreq(j_0_3*1e6)
    Omega_4 = Rabifreq(j_0_4*1e6)
    
    #Transition rates:
    dim = 7
    b = np.zeros(dim**2-1, dtype = complex)
    M = np.zeros([dim**2-1, dim**2-1], dtype = complex)
    Unit = np.zeros([dim**2-1, dim**2-1], dtype= complex) #Uninatry part of the master equation
    Diss = np.zeros([dim**2-1, dim**2-1], dtype= complex) #Dissipative part of the master equation
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
<<<<<<< HEAD
    Unit[dim*5+6, dim*5+5], Unit[dim*5+6, dim*4+6] = 1j*ks[5]*N(beta, ws[5])**0.5, -1j*ks[4]*N(beta, ws[4])**0.5
    Unit[dim*6+5, dim*6+4], Unit[dim*5+6, dim*5+5] = 1j*ks[4]*N(beta, ws[4])**0.5, -1j*ks[5]*N(beta, ws[5])**0.5
    for i in range(6):
        Unit[dim*5+6, dim*i+i] = +1j*ks[5]*N(beta, ws[5])**0.5
        Unit[dim*6+5, dim*i+i] = -1j*ks[5]*N(beta, ws[5])**0.5
    b[dim*5+6] += 1j*ks[5]*N(beta, ws[5])**0.5
    b[dim*6+5] -= 1j*ks[5]*N(beta, ws[5])**0.5
    Unit[dim*0+6, dim*0+5] = 1j*ks[5]*N(beta, ws[5])**0.5
    Unit[dim*0+6, dim*1+6] = -1j*ks[1]*N(beta, ws[1])**0.5
    Unit[dim*6+0, dim*6+1] = 1j*ks[1]*N(beta, ws[1])**0.5
    Unit[dim*6+0, dim*5+0] = -1j*ks[5]*N(beta, ws[5])**0.5
=======
        Unit[dim*5+6, dim*5+5], Unit[dim*5+6, dim*4+6] = 1j*ks[5]*N(beta, ws[5])**0.5, -1j*ks[4]*N(beta, ws[4])**0.5
        Unit[dim*6+5, dim*6+4], Unit[dim*5+6, dim*5+5] = 1j*ks[4]*N(beta, ws[4])**0.5, -1j*ks[5]*N(beta, ws[5])**0.5
        for i in range(6):
            Unit[dim*5+6, dim*i+i] = 1j*ks[5]*N(beta, ws[5])**0.5
            Unit[dim*5+6, dim*i+i] = 1j*ks[5]*N(beta, ws[5])**0.5
    
>>>>>>> abb4ce04aa2527b5933285903ea5a09e4494d735
    for i in range(dim-1):
        if i != 3:
            Diss[dim*i+i+1, dim*i+i+1] -= 0.5*gamma_nr
    for i in range(4):
        for j in range(4, dim):
            Diss[dim*i+j, dim*i+j] -= 0.5*gamma_nr
    for i in range(4):
        for j in range(4):
            Diss[dim*i+i, dim*j+j] -= gamma
<<<<<<< HEAD
    Diss[dim*4+4, dim*4+4] -= 4*gamma
    Diss[dim*5+5, dim*5+5] -= 4*gamma
    for i in range(dim-2):
        if i!=3:
            Diss[dim*i+i, dim*(i+1)+i+1] += gamma_nr
    for i in range(dim-1):
        Diss[dim*5+5, dim*i+i] -= gamma_nr
    b[dim*5+5] -= gamma_nr
=======
    Diss[dim*5+5, dim*5+5] -= 4*gamma
    Diss[dim*6+6, dim*6+6] -= 4*gamma
    for i in range(dim-1):
        if i!=3:
            Diss[dim*i+i, dim*(i+1)+i+1] += gamma_nr
    for i in range(dim-1):
        Diss[]
>>>>>>> abb4ce04aa2527b5933285903ea5a09e4494d735
    for i in range(1, dim-1):
        if 1!=4:
            Diss[dim*i+i, dim*i+i] -= gamma_nr
    
    for i in range(4):
<<<<<<< HEAD
            b[dim*i+i] -= gamma
      

    M = Unit + Diss
    
    rho = np.linalg.solve(M,b)
    diag_rho = [rho[i] for i in range(0,dim**2-1,8)]
    diag_rho.append(1-np.sum(diag_rho)) 
    print(f"rho: {diag_rho}")   
    

    return diag_rho
=======
            b[dim*i+i] = -gamma

    M = Unit[:48,:48] + Diss[:48,:48]

    rho = np.linalg.solve(M,b)
    
    print(np.diag(rho), 1-np.sum(np.diag(rho)))
    
    return 0 
>>>>>>> abb4ce04aa2527b5933285903ea5a09e4494d735

Steady_state(300,0,0.6)


