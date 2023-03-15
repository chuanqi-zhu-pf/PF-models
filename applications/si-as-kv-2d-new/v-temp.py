import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dfc = pd.read_csv(f"v-temp-1d.csv", header=None, delimiter=" ")
x = dfc.loc[0].values
# print(x)
y = dfc.loc[1].values
plt.scatter(x, y)
plt.grid()

# plt.xlabel("Position ($\\times$$10^{-10}$m)", fontsize=12)
# plt.ylabel("As Concentration (at.)", fontsize=12)
plt.show()
