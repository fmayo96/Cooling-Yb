import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters
plt.rcParams.update({'font.size': 7,
                     'lines.linewidth': 1,
                     'mathtext.default': 'regular',
                     'font.family': 'serif',
                     'figure.dpi': 300,
                     'grid.color': (0.5, 0.5, 0.5, 0.3)})
Parameters.d *= 2
Parameters.alpha_imp /= 6

cavity_list = ["Reasonable cavity R_90_ecs.npz",
               "Reasonable cavity R_95_ecs.npz",
               "Reasonable cavity R_96_ecs.npz",
               "Reasonable cavity R_97_ecs.npz",
               "Reasonable cavity R_98_ecs.npz",
               "Reasonable cavity R_99_ecs.npz",
               "Reasonable cavity R_995_ecs.npz",
               "Reasonable cavity R_999_ecs.npz",
               ]

cavity_list2 = ["Cavity b_0_ecs.npz",
                "Cavity b_1_ecs.npz",
                "Cavity b_2_ecs.npz",
                "Cavity b_3_ecs.npz",
                "Cavity b_4_ecs.npz",
               ]

N = 100
Ts = [300, 275, 250, 225, 200, 175, 150, 125, 100, 78]
j0 = np.linspace(0, 2, N)

alphas = [1/i for i in range(1,len(cavity_list)+1)]
alphas = list(reversed(alphas))

plt.figure(figsize=(4,3.7))
for c,cavity in enumerate(cavity_list):
  eff_maxpow = np.zeros(len(Ts))
  for n,T in enumerate(Ts):
    pow = np.zeros(N)
    pabs = np.zeros(N)
    for i in range(N):
      pow[i], pabs[i] = Power.Net_Power(T, j0[i], 0.0, cavity=cavity)
    idx = np.argmax(pow)
    eff_maxpow[n] = pow[idx] / pabs[idx]

  #np.save("mod_cav.npy", eff_maxpow)

  eff_nocav = np.load("eff_nocav.npy")
  #eff_mod = np.load("mod_cav.npy")

  plt.plot(Ts, eff_maxpow, ".-C0", alpha=alphas[c], label = r"$R=0.$" + (cavity.split(" ")[2]).split(".")[0].split("_")[1] + " Pumping: "+ r"$2\leftrightarrow 5$")
plt.plot(Ts, eff_nocav, "--k", label = "No cavity: Pumping " + r"$4\leftrightarrow 5$")
#plt.plot(Ts, eff_mod, ".-k", label = "Purcell No inhibition")
plt.legend(loc="lower right")
plt.xlabel("T (K)", fontsize=9)
plt.ylabel(r"$\eta$", fontsize=9)
plt.xticks(Ts)
plt.tight_layout()
plt.savefig("eff_vs_T_vs_R.pdf")
plt.show()