import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

N = 100
Ts = [300, 250, 200, 150]

j0 = np.linspace(0, 1.5, N)
Parameters.d2 = Parameters.d
COLORS = ['C0', 'C3', 'C4', 'C7']
plt.figure()
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  pow_purcell = np.zeros(N)
  for i in range(N):
    pow[i], _ = Power.Net_Power(T, 0.0, j0[i])
    pow_purcell[i], _ = Power.Net_Power(T, j0[i],0.0, purcell=True)
  print(f"T: {T} | maxPow: {np.max(pow)}")
  print(f"T: {T} | maxPow (Purcell): {np.max(pow_purcell)}")
  #np.save(f"pow_purcell54_secondlaser_T={T}", pow)
  plt.plot(j0, pow, '--',color=COLORS[n], alpha=0.5,linewidth=2)  
  plt.plot(j0, pow_purcell, color=COLORS[n], linewidth=2, label = f"T={T}")  
plt.ylabel(r"$P_{cool}$")
plt.xlabel(r"$j_0$", fontsize=12)
plt.ylim([-0.1, 400])
plt.legend()
plt.show()