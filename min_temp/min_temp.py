import numpy as np 
import matplotlib.pyplot as plt 

pi = np.pi
h = 6.62607015e-34  # in J s
hbar = h/(2*pi)
c = 2.99792458e10   # cm/s
kB = 1.380649e-23   # J/K
eps0 = 8.8541878128e-14    # F/cm = C/ (V cm)
# Dipolar moment
d = 3.31e-31 # C cm
#Transition energies in J (we converted the data from cm^-1 to J by multiplying by hc)
w = 9811*h*c
w1 = 237*h*c
w2 = 138*h*c
w3 = 102*h*c
w4 = 132*h*c
w5 = 150*h*c
#Spontaneous Emission in s^-1
gamma = 314.5
#Non resonant transition rates in s^-1
gamma_nr = 1e12

#Non radiative decay rate in s^-1
W_nr = 1.45 

#Refraction index
n = 1.45
#Absortion Coefficients
alpha_imp = 4e-4 # cm-1
alpha_rad = 1e-2 # cm-1
E1 = w1
E2 = E1 + w2
n_ion = 1.3e20 # cm-3
N_e = 13
eta_e = 0.92 # extraction efficiency


def N(T):
  return 1/(np.exp(w1/(kB*T)) - 1)

def Rabifreq(j_0):
    E_0 = np.sqrt(2*j_0*1e6/(c*n*eps0))
    return d*E_0/hbar

def Pnet(T, j_0):
  lambd = (gamma_nr * (N(T) + 1) - gamma * ((N(T)/(N(T)+1))**2 + (N(T)/(N(T)+1)) + 1)) / (gamma_nr * N(T) + gamma)
  rho_22 = ((2*hbar**2 * Rabifreq(j_0)**2)/(gamma_nr*(2*N(T)+1)) * (gamma_nr*N(T))/(gamma_nr*N(T)+gamma)
            + gamma*gamma_nr*N(T)/(gamma_nr*N(T)+gamma)) / (gamma_nr*(N(T)+1) + (2*hbar**2 * Rabifreq(j_0)**2)/(gamma_nr*(2*N(T)+1)) * (2+(N(T)/(N(T)+1))**2
                 +(N(T)/(N(T)+1))+lambd) - gamma_nr*lambd*N(T))
  rho_55 = rho_22 / (1 + gamma*gamma_nr*(2*N(T)+1) / (2*hbar**2*Rabifreq(j_0)**2))
  Pnet = gamma*rho_55*(E2-E1)* n_ion * N_e * eta_e - alpha_imp * j_0 * 1e6
  return Pnet 

Ts = [70,80,90,100]
Nj = 1000
j0 = np.linspace(0,1,Nj)
plt.figure()
for T in Ts:
  pow = np.zeros(Nj)
  for i in range(Nj):
     pow[i] = Pnet(T, j0[i])
  plt.plot(j0, pow, linewidth=2, label=r"$T=$" + str(T))
plt.legend()
plt.xlabel(r"$j_0$", fontsize=12)
plt.ylabel("Pow", fontsize=12)
plt.show()