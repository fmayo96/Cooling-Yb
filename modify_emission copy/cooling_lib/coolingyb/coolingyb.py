from qutip import *
import numpy as np
from .aux_func import *
from .steady_state import Steady_state
from .parameters import Parameters as P


def Lindbladian(c_op : Qobj, rho : Qobj) -> Qobj:
  L = c_op * rho * c_op.dag() - 0.5 * (c_op.dag() * c_op * rho + rho * c_op.dag() * c_op )
  return L

class Power:
  def Net_Power(T, j_0_3, j_0_4, purcell=False,heating_data=False):
    rho = Steady_state(T, j_0_3, j_0_4, purcell)
    
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

    c_decay = np.zeros((7,7))
    for i in range(4):
        for j in range(4, 7):
            c_decay[i,j] = 1
    c_decay = np.sqrt(P.W_nr) * c_decay
    c_decay = Qobj(c_decay)
    
    L = Lindbladian(c_nr, rho) 
    L +=  Lindbladian(c_nr_dag, rho) 
    L += Lindbladian(c_decay, rho)

    H_laser = Qobj(Hamiltonian(j_0_3, j_0_4))
    unitary = -1j*(H_laser * rho - rho * H_laser)
    n_ion = 1.3e20 # cm-3
    N_e = 13
    eta_e = 0.92 # extraction efficiency

    Hs = np.diag([P.E0, P.E1, P.E2, P.E3, P.E4, P.E5, P.E6])
    Hs = Qobj(Hs/P.hbar)
    Pimp = P.alpha_imp * (j_0_4 + j_0_3) * 1e6
    Pcool = (Hs * L).tr() * n_ion * N_e * eta_e * P.hbar - Pimp
    #Pabs = (j_0_4 + j_0_3) * 1e6 * (P.alpha_imp + P.alpha_rad)
    Pabs_me = (Hs * unitary).tr() * n_ion * N_e * eta_e * P.hbar + Pimp

    if heating_data == True:
      non_rad_heating = -(Hs * Lindbladian(c_decay, rho)).tr() * n_ion * N_e * eta_e * P.hbar
      return np.real(Pcool), np.real(Pabs_me), np.real(non_rad_heating), np.real(Pimp)   

    return np.real(Pcool), np.real(Pabs_me)
