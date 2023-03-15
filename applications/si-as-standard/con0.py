import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.figure(figsize=(10, 8))

fig, ax = plt.subplots()

nx = 128
ny = 128
step = 12000
x = np.linspace(0, 10*ny, num=ny)
dfc = pd.read_csv(f"128-2d/conl/2d{step}.csv", header=None)
arrc = dfc[0].values
matc = arrc.reshape(ny, nx)
ax.plot(x, matc[:, 0], linewidth=2.0, label="$c_L^{\eta_2=70d_0}$")

nx = 256
ny = 256
step = 48000
x = np.linspace(0, 5*ny, num=ny)
dfc = pd.read_csv(f"256-2d/conl/2d{step}.csv", header=None)
arrc = dfc[0].values
matc = arrc.reshape(ny, nx)
ax.plot(x, matc[:, 0], linewidth=1.0, label="$c_L^{\eta_3=35d_0}$")

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
# ax.legend(fontsize=18, bbox_to_anchor=(1, 0.9))

ax2 = ax.twinx()
nx = 128
ny = 128
step = 12000
x = np.linspace(0, 10*ny, num=ny)
dfc = pd.read_csv(f"128-2d/phi/2d{step}.csv", header=None)
arrc = dfc[0].values
matc = arrc.reshape(ny, nx)
ax2.plot(x, matc[:, 0], linestyle="dashed",
         linewidth=2.0, label="$\phi_L^{\eta_2=70d_0}$")

nx = 256
ny = 256
step = 48000
x = np.linspace(0, 5*ny, num=ny)
dfc = pd.read_csv(f"256-2d/phi/2d{step}.csv", header=None)
arrc = dfc[0].values
matc = arrc.reshape(ny, nx)
ax2.plot(x, matc[:, 0],  linestyle="dashed",
         linewidth=1.0, label="$\phi_L^{\eta_3=35d_0}$")

# ax2.legend(fontsize=18, bbox_to_anchor=(1, 0.7))

# plt.xlim([410, 470])
# plt.ylim([0.0155, 0.0182])

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.show()
