import numpy as np 
import matplotlib.pyplot as plt 

with np.load("Intermediate cavity for 971 nm_abs.npz",allow_pickle=True) as f:
    convo_lines_T = f["convo_lines_T"]
    lambda_lines = f["lambda_lines"]

plt.figure()
plt.plot(lambda_lines, convo_lines_T[-1], '.')
plt.xlabel("wavelength", fontsize=12)
plt.show()
