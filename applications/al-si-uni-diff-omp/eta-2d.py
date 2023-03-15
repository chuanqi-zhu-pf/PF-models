import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 2000
nx = 128
ny = 128
step_arr = np.arange(0, ns*11, ns)
for step in step_arr:
    dfc = pd.read_csv(f"128dx/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.plot(matc[0, :])
plt.show()

# dfc = pd.read_csv(f"128-conp0-5/con/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.plot(matc[0, :])
