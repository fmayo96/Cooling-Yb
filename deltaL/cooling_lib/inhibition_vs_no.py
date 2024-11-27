import numpy as np 
import matplotlib.pyplot as plt 
from coolingyb import Power, Parameters

Parameters.d *= 2
Parameters.alpha_imp /= 6


CAVITY = "Reasonable cavity R_98_ecs.npz"
CAVITY_NOINHIB = "Reasonable cavity R_98_ecs_noinhib.npz"

N = 100
Ts = [300, 275, 250, 225, 200, 175, 150, 125, 100, 78]
j0 = np.linspace(0, 2, N)

pow = np.zeros(N)
pabs = np.zeros(N)
pow_no = np.zeros(N)
pabs_no = np.zeros(N)
eff = np.zeros(N)
eff_no = np.zeros(N)
effMax = np.zeros(len(Ts))
effMax_no = np.zeros(len(Ts))
for n, T in enumerate(Ts):
  for i in range(N):
    pow[i], pabs[i] = Power.Net_Power(T, j0[i], 0, cavity=CAVITY)
    eff[i] = pow[i] / pabs[i]
    pow_no[i], pabs_no[i] = Power.Net_Power(T, j0[i], 0, cavity=CAVITY_NOINHIB)
    eff_no[i] = pow_no[i] / pabs_no[i]
  idx = np.argmax(pow)
  effMax[n] = eff[idx]
  idx = np.argmax(pow_no)
  effMax_no[n] = eff_no[idx]


plt.figure()
plt.plot(Ts, effMax/effMax_no, '.-', label = "Inhibition")
#plt.plot(Ts, effMax_no, '.-k', label = "No Inhibition")
plt.xlabel('T (K)', fontsize=12)
plt.ylabel(r"$\eta$", fontsize=12)
plt.legend(fontsize=11)
#plt.savefig("Inhibition_vs_NoInhibition.pdf", dpi=800)
plt.show()