import numpy as np 
import matplotlib.pyplot as plt 

with np.load("lines.npz") as data:
  lambda_lines = data["lambda_lines"]
  ECS_E_c_lines_T = data["ECS_E_c_lines_T"]
  ACS_E_c_lines_T = data["ACS_E_c_lines_T"]



plt.figure()
plt.plot(lambda_lines, ECS_E_c_lines_T[-1], '.', label = "Emision")
plt.plot(lambda_lines, ACS_E_c_lines_T[-1], '.', label = "Absortion")
plt.legend()
plt.xlabel(r"$\lambda$", fontsize=12)
plt.show()

