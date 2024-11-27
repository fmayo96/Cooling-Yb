import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters


cavity_list = ["Feld 1987 for 971 nm.npz",
               "Intermediate cavity for 971 nm.npz", 
               "Spherical cavity for 971 nm.npz"]


N = 100
Ts = [300, 250, 200, 150, 125, 100, 78]
j0 = np.linspace(0, 2, N)

COLORS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5','C6'] #['C3', 'C4', 'C0', 'C7', 'C8', 'C9','k']

plt.figure()
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  pabs = np.zeros(N)
  for i in range(N):
    pow[i], pabs[i] = Power.Net_Power(T, j0[i], 0.0, cavity=cavity_list[1])
  plt.plot(j0, pow/pabs, color=COLORS[n], linewidth=2, label = f"T={T}")  
plt.ylabel(r"$\eta$", fontsize=12)
plt.xlabel(r"$j_0\,\,(MW/cm^2)$", fontsize=12)
#plt.ylim([-0.1, 500])
plt.legend()
plt.savefig("efficiency_feld_cav", dpi=800)
plt.show()

