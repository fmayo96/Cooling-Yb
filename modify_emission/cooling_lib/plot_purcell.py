import numpy as np 
import matplotlib.pyplot as plt

Ts = [300, 275, 250, 225, 200, 175, 150]
N = 100

pow_nopurcell = np.zeros((len(Ts), N))
for i,T in enumerate(Ts):
  pow_nopurcell[i] = np.load(f'pow_nopurcell_T={T}.npy')

pow_purcell54 = np.zeros((len(Ts), N))
for i,T in enumerate(Ts):
  pow_purcell54[i] = np.load(f'pow_purcell54_T={T}.npy')

pow_purcell54_secondlaser = np.zeros((len(Ts), N))
for i,T in enumerate(Ts):
  pow_purcell54_secondlaser[i] = np.load(f'pow_purcell54_secondlaser_T={T}.npy')

pow_max_nopurcell = np.array([np.max(pow_nopurcell[i]) for i in range(len(Ts))])
pow_max_purcell = np.array([np.max(pow_purcell54[i]) for i in range(len(Ts))])
pow_max_purcell_secondlaser = np.array([np.max(pow_purcell54_secondlaser[i]) for i in range(len(Ts))])



plt.figure()
plt.plot(Ts, pow_max_purcell/pow_max_nopurcell, '.', label="laser 45")
plt.plot(Ts, pow_max_purcell_secondlaser/pow_max_nopurcell, '.', label='laser 35')
plt.plot(Ts, np.ones(len(Ts)), '--k')
plt.xlabel("Temperature", fontsize=12)
plt.ylabel(r"$P_{max}^{(purcell)} / P_{max}$", fontsize=12)
plt.legend()
plt.show()