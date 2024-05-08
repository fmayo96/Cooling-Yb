from qutip import *
import numpy as np
from parameters import *
from aux_func import *
from steady_state import Steady_state
import matplotlib.pyplot as plt 


def Lindbladian(c_op : Qobj, rho : Qobj) -> Qobj:
  L = c_op * rho * c_op.dag() - 0.5 * (c_op.dag() * c_op * rho + rho * c_op.dag() * c_op * rho)
  return L



def Net_Power(T, j_0_3, j_0_4):
  rho = Steady_state(T, j_0_3, j_0_4)
  
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

  #L = Lindbladian(c_nr, rho)
  #L += Lindbladian(c_nr.dag(), rho)

  rho_ss = Steady_state(300, 0.0, 0.0)
  L = Lindbladian(c_nr, rho_ss)
  L += Lindbladian(c_nr.dag(), rho_ss)
  L += Lindbladian(c_se, rho_ss)
  print(L)


  n_ion = 1.3e20 # cm-3
  N_e = 13
  eta_e = 0.92 # extraction efficiency

  Hs = np.diag([E0, E1, E2, E3, E4, E5, E6])
  Hs = Qobj(Hs)
  Pcool = (Hs * L).tr() * n_ion * N_e * eta_e / 1e10  - alpha_imp * j_0_4 * 1e6
  return Pcool

Net_Power(300, 0, 0)
# Ts = [300, 150, 110, 100]
# Ns = 1000
# j_0_4 = np.linspace(0, 1.8, Ns)
# plt.figure()
# for T in Ts:
#    pow = np.zeros(Ns)
#    for i in range(Ns):
#       pow[i] = Net_Power(T, 0.0, j_0_4[i])
#    pow = np.array(pow)
#    #plt.subplot(1,2,1)
#    plt.plot(j_0_4, pow-pow[0], linewidth = 2, label = r'$T=$' + str(T))
#    #plt.subplot(1,2,2)
#    #plt.plot(j_0_4, eta, linewidth = 2, label = r'$T=$' + str(T))
# plt.legend()
# plt.ylim([0, 300])
# plt.show()



