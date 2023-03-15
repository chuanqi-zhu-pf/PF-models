import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.figure(figsize=(12, 9), dpi=80)
dfc = pd.read_csv(f"kv-1603-8dx/con/1d200000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200])
dfc = pd.read_csv(f"kv-1580-8dx/con/1d200000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200])
dfc = pd.read_csv(f"kv-1560-8dx/con/1d200000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200])
plt.xlabel("Position ($\\times$$10^{-10}$m)", fontsize=12)
plt.ylabel("As Concentration (at.)", fontsize=12)
plt.show()
