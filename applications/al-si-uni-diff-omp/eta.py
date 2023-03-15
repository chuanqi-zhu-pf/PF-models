import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

step = 1200
dfc = pd.read_csv(f"5dx/con/1d{step}.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc)

dfc = pd.read_csv(f"10dx/con/1d{step}.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc)
plt.show()
