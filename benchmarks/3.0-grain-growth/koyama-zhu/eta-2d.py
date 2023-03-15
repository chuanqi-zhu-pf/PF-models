import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 128
ny = 128
step = 1000
dfc = pd.read_csv(f"7dx/2d{step}.csv", header=None)
arrc = dfc[0].values
matc = arrc.reshape(ny, nx)
plt.plot(matc[:, 63])


dfc = pd.read_csv(f"4dx/2d{step}.csv", header=None)
arrc = dfc[0].values
matc = arrc.reshape(ny, nx)
plt.plot(matc[:, 63])
plt.show()
