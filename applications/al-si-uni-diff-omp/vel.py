import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import scipy as sc

# dfv = pd.read_csv(f"data/vel.csv", header=None)
# dft = pd.read_csv(f"data/time.csv", header=None)
# x = dft[0].values
# y = dfv[0].values

# fit fourth-degree polynomial
# model4 = np.poly1d(np.polyfit(x, y, 18))

# define scatterplot
# polyline = np.linspace(0.001494, 0.29, 1000)
# plt.plot(x, y)

dfv = pd.read_csv(f"64-5dx/vel.csv", header=None)
dft = pd.read_csv(f"64-5dx/time.csv", header=None)
x = dft[0].values
y = dfv[0].values

plt.plot(x, y)

# add fitted polynomial curve to scatterplot
# plt.plot(polyline, model4(polyline), "--", color='red', linewidth=1)
plt.show()

# plt.show()

# plt.plot(time, vel)

# plt.show()
