import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

Nstep = 1000
Ts = [150, 140, 138, 130]
js = np.linspace(0, 0.2, Nstep)
COLORS = ['C0', 'C1', 'C2', 'C3']

plt.figure()
for n,T in enumerate(Ts):
  pow = np.zeros(Nstep)  
  pow_sl = np.zeros(Nstep)
  for i,j0 in enumerate(js):
    pow[i], _, _ = Power.Net_Power(T, 0, j0)
    pow_sl[i], _, _ = Power.Net_Power(T, j0, 0)
  plt.plot(js, pow, color=COLORS[n], linewidth = 2, label = r"$T=$" + str(T))
  plt.plot(js, pow_sl, '--', color=COLORS[n], linewidth = 2)
plt.legend(fontsize=11)
plt.xlabel(r"$\mathrm{Pump\,\,\,intensity}\,MW/cm^2$", fontsize=12)
plt.ylabel(r"$\mathrm{Net\,\,\,Cooling\,\,\,Power}\,W/cm^3$", fontsize=12)
plt.plot(js, np.zeros(Nstep), '--k')
plt.ylim([-0.1, 10])
plt.show()