import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power


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

plt.figure()
for cavity in cavity_list:
  eff_maxpow = np.zeros(len(Ts))
  for n,T in enumerate(Ts):
    pow = np.zeros(N)
    pabs = np.zeros(N)
    for i in range(N):
      pow[i], pabs[i] = Power.Net_Power(T, j0[i], 0.0, cavity=cavity)
    idx = np.argmax(pow)
    eff_maxpow[n] = pow[idx] / pabs[idx]

  #np.save("mod_cav.npy", eff_maxpow)

  eff_nocav = np.load("eff_maxpow.npy")
  #eff_mod = np.load("mod_cav.npy")

  plt.plot(Ts, eff_maxpow, ".-", label = (cavity.split(" ")[2]).split(".")[0])
plt.plot(Ts, eff_nocav, "--k", label = "No cavity")
#plt.plot(Ts, eff_mod, ".-k", label = "Purcell No inhibition")
plt.legend()
plt.xlabel("T (K)", fontsize=14)
plt.ylabel(r"$\eta$", fontsize=14)
plt.xticks(Ts, fontsize=11)
plt.yticks(fontsize=11)
plt.savefig("eff_vs_T.png", dpi=800)
plt.show()