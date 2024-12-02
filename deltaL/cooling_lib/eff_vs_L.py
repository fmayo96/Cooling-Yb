import numpy as np 
from coolingyb import Power, Parameters
import matplotlib.pyplot as plt 
plt.rcParams.update({'font.size': 7,
                     'lines.linewidth': 1,
                     'mathtext.default': 'regular',
                     'font.family': 'serif',
                     'figure.dpi': 300,
                     'grid.color': (0.5, 0.5, 0.5, 0.3)})

Parameters.d *= 2
Parameters.alpha_imp /= 6

Ls = np.array([79965 + i*5 for i in range(18)])

cavity_list_96 = [f"Reasonable cavity R_96 L_{L}_ecs.npz" for L in Ls]
cavity_list_99 = [f"Reasonable cavity R_99 L_{L}_ecs.npz" for L in Ls]
cavity_list = [cavity_list_96, cavity_list_99]
R = [0.96, 0.99]
COLORS = ['C0', 'k']
T = 100
N = 100
j0 = np.linspace(0, 2, N)

plt.figure(figsize=(3.5, 3.3))
for n, cavities in enumerate(cavity_list):
  eff_max_pow = []
  for cavity in cavities:
    pow, pabs = np.zeros(N), np.zeros(N)
    for i in range(N):
      pow[i], pabs[i] = Power.Net_Power(T, j0[i], 0.0, cavity=cavity)
      idx = np.argmax(pow)
    eff_max_pow.append(pow[idx]/pabs[idx])
  plt.plot((Ls[:-2] - 80000)/10, eff_max_pow[:-2], '.-', color=COLORS[n], label = r"$R=$" + str(R[n]))
plt.legend()
plt.xlabel(r"$\delta L$" + " (nm)", fontsize=9)
plt.ylabel(r"$\eta$", fontsize=9)
plt.tight_layout()
plt.xticks([-4,-3,-2,-1,0,1,2,3,4])
plt.savefig("eff_vs_deltaL.pdf")
plt.show()