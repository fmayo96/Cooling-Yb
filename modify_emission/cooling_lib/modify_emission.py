import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

N = 300
Ts = [300, 275, 250, 225, 200, 175, 150, 125, 100, 78]
j0 = np.linspace(0, 2, N)


for n,T in enumerate(Ts):
  pow = np.zeros(N)
  pabs = np.zeros(N)
  eff_maxpow = np.zeros(len(Ts))
  for i in range(N):
    pow[i], pabs[i] = Power.Net_Power(T, 0.0, j0[i])
  plt.plot(j0, pow/pabs,  linewidth=2, label = f"T={T}")  
plt.ylabel(r"$\eta$", fontsize=12)
plt.xlabel(r"$j_0\,\,(MW/cm^2)$", fontsize=12)
plt.ylim([-0.001, 0.015])
plt.legend()
plt.savefig("eff_no_cavity", dpi=800)
plt.show()
