import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.figure(figsize=(12, 9), dpi=80)
dfc = pd.read_csv(f"kv-1603-8dx/con/1d200000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black", label="1603 K",linewidth=0.5)
dfc = pd.read_csv(f"kv-1580-8dx/con/1d200000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200],color="black", label="1580 K", linewidth=1.0)
dfc = pd.read_csv(f"kv-1560-8dx/con/1d200000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black", label="1560 K",  linewidth=1.5)
plt.xlabel("Position ($\\times$$10^{-10}$m)", fontsize=14)
plt.ylabel("As Concentration (at.)", fontsize=14)
plt.legend(fontsize=14)
plt.show()
