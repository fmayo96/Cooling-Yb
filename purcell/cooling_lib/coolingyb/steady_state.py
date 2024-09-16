import numpy as np
from qutip import *
from .parameters import Parameters as P
from .aux_func import *


ORDER = [4, 7, 9, 11, 1, 5, 6, 10, 0, 2, 3, 8]

with np.load('Purcell_enhanced_emission.npz') as data:
    g = data['convo_lines_T_n']
    g = g[1][3]
    g = g[:,ORDER]
    g = g.reshape(10, 3, 4)
Ts = [78, 100, 125, 150, 175, 200, 225, 250, 275, 300]
gammas = {}
for i,T in enumerate(Ts):
    gammas[T] = g[i]



def Steady_state(T, j_0_3, j_0_4, inhib):
    beta = 1/(T*P.kB)
    c_nr = np.zeros((7,7))
    ws = [P.w1, P.w2, P.w3, P.w4, P.w4, P.w5]
    for i in range(6):
        c_nr[i,i+1] = np.sqrt(N_therm(beta, ws[i]) + 1)
    c_nr[3,4] = 0
    c_nr = np.sqrt(P.gamma_nr)* c_nr

    c_nr = Qobj(c_nr)

    c_nr_dag = np.zeros((7,7))
    for i in range(1, 7):
        c_nr_dag[i,i-1] = np.sqrt(N_therm(beta, ws[i-1]))
    c_nr_dag[4,3] = 0

    c_nr_dag = Qobj(c_nr_dag * np.sqrt(P.gamma_nr))

    c_se = np.zeros((7,7))
    for i in range(4):
        for j in range(4, 7):
            c_se[i,j] = (max(gammas[T][j-4,i], 0.0))**0.5
            #c_se[i,j] = 1
    
    #c_se[3,4] *= (1-inhib) 
    #c_se[2,4] *= (1-inhib)
    #c_se[1,4] = 0
    c_se *= P.gamma**0.5*1.5
    c_se = Qobj(c_se)

    c_decay = np.zeros((7,7))
    for i in range(4):
        for j in range(4, 7):
            c_decay[i,j] = 1
    c_decay = np.sqrt(P.W_nr) * c_decay
    c_decay = Qobj(c_decay)

    Hs = np.diag([P.E0, P.E1, P.E2, P.E3, P.E4, P.E5, P.E6])
    Hs = Qobj(Hs)

    H = Hamiltonian(j_0_3, j_0_4)
    H = Qobj(H) 

    Htot = H #+ Hs
    rho_ss = steadystate(Htot, [c_nr, c_nr_dag, c_se, c_decay])
    return rho_ss


