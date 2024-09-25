import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

N = 300
Ts = [300, 250, 200, 150, 100, 78]
j0 = np.linspace(0, 10, N)

COLORS = ['C3', 'C4', 'C0', 'C7', 'C8', 'k']

plt.figure()
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  for i in range(N):
    pow[i], _ = Power.Net_Power(T, 0.0, j0[i])
  plt.plot(j0, pow, color=COLORS[n], linewidth=2, label = f"T={T}")  
plt.ylabel(r"$P_{cool} (MW/cm^3)$", fontsize=12)
plt.xlabel(r"$j_0\,\,(MW/cm^2)$", fontsize=12)
plt.ylim([-0.1, 500])
plt.legend()
plt.savefig("P_cool", dpi=800)
plt.show()
