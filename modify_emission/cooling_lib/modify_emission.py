import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

N = 100
Ts = [150, 125, 100]
j0 = np.linspace(0, 1.5, N)
inhib = np.linspace(0,1,N)

COLORS = ['C3', 'C4', 'C0', 'C7']

plt.figure()
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  pow_purcell = np.zeros(N)
  max_pow = np.zeros(N)
  max_pow_purcell = np.zeros(N)
  for j in range(N):
    for i in range(N):
      pow[i], _ = Power.Net_Power(T, 0.0, j0[i])
      pow_purcell[i], _ = Power.Net_Power(T, j0[i],0.0, inhib[j])
      max_pow[j] = max(pow)
      max_pow_purcell[j] = max(pow_purcell)
  plt.plot(inhib, max_pow, '--',color=COLORS[n], alpha=0.5,linewidth=2)  
  plt.plot(inhib, max_pow_purcell, color=COLORS[n], linewidth=2, label = f"T={T}")  
plt.ylabel(r"$P_{cool} (MW/cm^3)$", fontsize=12)
plt.xlabel(r"$\mathrm{Purcell}\,\,\,\mathrm{Inhibition}$", fontsize=12)
#plt.ylim([-0.1, 400])
plt.legend()
#plt.savefig("P_cool_inhibition_25", dpi=800)
plt.show()