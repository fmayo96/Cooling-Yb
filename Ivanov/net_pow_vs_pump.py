import numpy as np
import matplotlib.pyplot as plt 
from power import Net_Power 

Temps = [300,150,100,70]

j_0_4 = np.linspace(0,1.8,20)

P_net = np.zeros([len(Temps), len(j_0_4)])

for i in range(len(Temps)):
    for j in range(len(j_0_4)):
        P_net[i,j] = Net_Power(Temps[i], 0, j_0_4[j])


plt.figure()
plt.plot(j_0_4, -P_net[0,:], linewidth = 2, label = 'T=' + str(Temps[0]))
plt.plot(j_0_4, -P_net[1,:], linewidth = 2, label = 'T=' + str(Temps[1]))
plt.plot(j_0_4, -P_net[2,:], linewidth = 2, label = 'T=' + str(Temps[2]))
plt.plot(j_0_4, -P_net[3,:], linewidth = 2, label = 'T=' + str(Temps[3]))
plt.legend()
plt.xlabel('Pump intensity ' + r'$MW/cm^2$', fontsize = 14)
plt.ylabel('Net Cooling Power ' + r'$MW/cm^3$', fontsize = 14)
#plt.ylim([0,400])
plt.show()  