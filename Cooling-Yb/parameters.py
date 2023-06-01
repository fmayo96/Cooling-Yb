#This File contains all the parameters used in the program.
#Physical constants, transition rates, etc.

import numpy as np 

# Universal constants
pi = np.pi
h = 6.62607015e-34  # in J s
hbar = h/(2*pi)
c = 2.99792458e10   # cm/s
kB = 1.380649e-23   # J/K
eps0 = 8.854e-14    # F/m = C/ (V cm)
# Dipolar moment
d = 3.34e-31 # C cm

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
#Refraction index
n = 1.45


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

H = np.array([E0,E1,E2,E3,E4,E5,E6])


