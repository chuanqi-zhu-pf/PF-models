import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dfc = pd.read_csv(f"tip_vel-16.csv", header=None)
arrc = dfc[0].values

Nx = len(arrc)
temp = dfc[0].values
temp2 = np.zeros(Nx)
dtime = 0.1
nstep = 200

for istep in range(nstep):
    for i in range(1, Nx-1):
        temp2[i] = temp[i] + dtime*(temp[i+1] + temp[i-1] - 2.0*temp[i])
    for i in range(1, Nx-1):
        temp[i] = temp2[i]

plt.plot(temp, label="16dx")

dfc2 = pd.read_csv(f"tip_vel-8.csv", header=None)
arrc2 = dfc2[0].values

Nx = len(arrc2)
temp = dfc2[0].values
temp2 = np.zeros(Nx)
dtime = 0.1
nstep = 200

for istep in range(nstep):
    for i in range(1, Nx-1):
        temp2[i] = temp[i] + dtime*(temp[i+1] + temp[i-1] - 2.0*temp[i])
    for i in range(1, Nx-1):
        temp[i] = temp2[i]

plt.plot(temp, label="8dx")
plt.legend()

plt.show()
