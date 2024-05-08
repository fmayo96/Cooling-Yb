#This File contains all the parameters used in the program.
#Physical constants, transition rates, etc.

import numpy as np 

# Universal constants
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
#Vibronic couplings in s^-1 (k/hbar)
k_12 = 6.35e11
k_23 = 2.8e11 
k_34 = 2.25e11
k_56 = 5.4e11 
k_67 = 4.1e11 

#Refraction index
n = 1.45
#Absortion Coefficients
alpha_imp = 4e-4 # cm-1
alpha_rad = 1e-2 # cm-1
#Energy levels
E0 = 0
E1 = w1
E2 = E1 + w2
E3 = E2 + w3
E4 = E3 + w
E5 = E4 + w4
E6 = E5 + w5 
E_gs = np.array([E0,E1,E2,E3])
E_es = np.array([E4,E5,E6])

H = np.zeros(7, dtype=np.complex128)
H[0] = E0
H[1] = E1
H[2] = E2
H[3] = E3
H[4] = E4
H[5] = E5
H[6] = E6


# print("w1 = ", w1)
# print("k1 = ", k_12 * hbar)
# print(k_12 ** 2 * hbar)
