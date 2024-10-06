import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power


cavity_list = ["Feld 1987 for 971 nm.npz",
               "Intermediate cavity for 971 nm.npz", 
               "Spherical cavity for 971 nm.npz"]


N = 100
Ts = [300, 275, 250, 225, 200, 175, 150, 125, 100, 78]
j0 = np.linspace(0, 2, N)

eff_maxpow = np.zeros(len(Ts))
for n,T in enumerate(Ts):
  pow = np.zeros(N)
  pabs = np.zeros(N)
  for i in range(N):
    pow[i], pabs[i] = Power.Net_Power(T, j0[i], 0.0, cavity=cavity_list[1])
  idx = np.argmax(pow)
  eff_maxpow[n] = pow[idx] / pabs[idx]

#np.save("mod_cav.npy", eff_maxpow)

eff_nocav = np.load("eff_maxpow.npy")
eff_mod = np.load("mod_cav.npy")

plt.figure()
plt.plot(Ts, eff_nocav, ".-C0", label = "No cavity")
plt.plot(Ts, eff_maxpow, ".-C3", label = "Purcell (La nuestra)")
plt.plot(Ts, eff_mod, ".-k", label = "Purcell No inhibition")
plt.legend()
plt.xlabel("T (K)", fontsize=14)
plt.ylabel(r"$\eta$", fontsize=14)
plt.xticks(Ts, fontsize=11)
plt.yticks(fontsize=11)
plt.savefig("eff_vs_T.png", dpi=800)
plt.show()