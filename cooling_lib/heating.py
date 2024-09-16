import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power

N = 100
Ts = [300, 250, 200, 150]
COLORS = ["C0", "C1", "C2", "C3"]
js = np.linspace(0, 1, N)

plt.figure()
for i, T in enumerate(Ts):
  pow = np.zeros(N)
  pabs = np.zeros(N)
  nr_heating = np.zeros(N)
  pimp = np.zeros(N)
  for j in range(N):
    pow[j], pabs[j], nr_heating[j], pimp[j] = Power.Net_Power(T, 0, js[j], heating_data=True)
  plt.plot(js, nr_heating, linewidth = 2, label = "T=" + str(T))
plt.plot(js, pimp, '--k', linewidth = 2, label = r"$P_{imp}$")
plt.xlabel(r"$\mathrm{Pump\,\,\,Intensity}\,\,MW/cm^2$", fontsize=12)
plt.ylabel(r"$\dot{Q}_{nr}$", fontsize = 12)
plt.legend()
plt.savefig("../results/nr_heating/nr_vs_pimp.png", dpi=300)
plt.show()