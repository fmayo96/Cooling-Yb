from qutip import *
from parameters import *
from aux_func import *

def Steady_state(T, j_0_3, j_0_4):
    beta = 1/(T*kB)
    c_nr = np.zeros((7,7))
    ws = [w1, w2, w3, w4, w4, w5]
    for i in range(6):
        c_nr[i,i+1] = np.sqrt(N(beta, ws[i]) + 1)
    c_nr[3,4] = 0
    c_nr = np.sqrt(gamma_nr/10)* c_nr

    c_nr = Qobj(c_nr)

    c_nr_dag = np.zeros((7,7))
    for i in range(1, 7):
        c_nr_dag[i,i-1] = np.sqrt(N(beta, ws[i-1]))
    c_nr_dag[4,3] = 0

    c_nr_dag = Qobj(c_nr_dag * np.sqrt(gamma_nr/10))

    c_se = np.zeros((7,7))
    for i in range(4):
        for j in range(4, 7):
            c_se[i,j] = 1
    c_se = np.sqrt(gamma)*c_se
    c_se = Qobj(c_se)


    H = Hamiltonian(j_0_3, j_0_4)
    H = Qobj(H) 
    rho_ss = steadystate(H, [c_nr, c_nr_dag, c_se])
    return rho_ss


