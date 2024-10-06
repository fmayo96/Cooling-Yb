import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

N = 300
Ts = [300, 275, 250, 225, 200, 175, 150, 125, 100, 78]
j0 = np.linspace(0, 2, N)

eff_maxpow = np.zeros(len(Ts))
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  pabs = np.zeros(N)
  for i in range(N):
    pow[i], pabs[i] = Power.Net_Power(T, 0, j0[i])
  idx = np.argmax(pow)
  eff_maxpow[n] = pow[idx] / pabs[idx]

np.save("eff_maxpow",eff_maxpow)

plt.figure()
plt.plot(Ts, eff_maxpow, '.')
plt.xlabel("T (K)", fontsize=14)
plt.ylabel(r"$\eta$", fontsize=14)
plt.show()