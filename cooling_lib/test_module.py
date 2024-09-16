import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters
from matplotlib import cm
from matplotlib import colors

N = 100
Ts = [300, 250, 200, 150]

j34 = np.linspace(0, 1.5, N)
j24 = np.linspace(0, 1.5, N)

J24, J34 = np.meshgrid(j34, j34)

plt.figure(figsize=(10,10))
for n,T in enumerate(Ts):
  Pow = np.zeros((N, N))
  Pabs = np.zeros((N, N))
  eff = np.zeros((N, N))

  for i in range(N):
    for j in range(N):
      Pow[i,j], Pabs[i,j] = Power.Net_Power(T, j24[i], j34[j])

  for i in range(N):
    for j in range(N):
      if Pabs[i,j] == 0:
        eff[i,j] = 0
      else:
        eff[i,j] = Pow[i,j] / Pabs[i,j]

  plt.subplot(2,2, n+1)
  plt.contourf(J34, J24, eff.T, 1000, cmap=cm.RdBu, vmax = np.max(eff), vmin = -np.max(eff))
  plt.title("Cooling efficiency T = " + str(T))
  if n >= 2:
    plt.xlabel(r"$\mathrm{Pump\,\,\,Intensity\,\,\,}4\leftrightarrow 5\,\,(MW/cm^2)$", fontsize=12)
  if n %2 ==  0:
    plt.ylabel(r"$\mathrm{Pump\,\,\,Intensity\,\,\,}3\leftrightarrow 5\,\,(MW/cm^2)$", fontsize=12)
  plt.colorbar()
#plt.savefig('../results/two_laser_plots/eff_subplots.png', dpi=300)
plt.show()


# plt.figure()
# plt.contourf(J34, J24, Pow.T, 1000, cmap=cm.RdBu, vmax = np.max(Pow), vmin = -np.max(Pow))
# plt.title(r"$\mathrm{Net\,\,\,Cooling\,\,\,Power}\,\,\,W/cm^3\,\,\,T = $" + str(T) + r"$\,K$")
# plt.xlabel(r"$\mathrm{Pump\,\,\,Intensity\,\,\,}4\leftrightarrow 5\,\,(MW/cm^2)$", fontsize=12)
# plt.ylabel(r"$\mathrm{Pump\,\,\,Intensity\,\,\,}3\leftrightarrow 5\,\,(MW/cm^2)$", fontsize=12)
# plt.colorbar()
# plt.show()

