import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

N = 200
Ts = [300, 250, 200, 150, 100, 78]
j0 = np.linspace(0, 5, N)

COLORS = ['C3', 'C4', 'C0', 'C7', 'C8', 'k']

plt.figure()
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  for i in range(N):
    pow[i], _ = Power.Net_Power(T, j0[i], 0.0)
  plt.plot(j0, pow, color=COLORS[n], linewidth=2, label = f"T={T}")  
plt.ylabel(r"$P_{cool} (MW/cm^3)$", fontsize=12)
plt.xlabel(r"$j_0$", fontsize=12)
plt.ylim([-0.1, 500])
plt.legend()
#plt.savefig("P_cool_inhibition_25", dpi=800)
plt.show()

