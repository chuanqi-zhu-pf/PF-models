import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"64-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]

plt.plot(x, y, label = "w/$d_0$=50")

df = pd.read_csv(f"128-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]

plt.plot(x, y, label = "w/$d_0$=25")

df = pd.read_csv(f"tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]

plt.plot(x, y, label = "w/$d_0$=12.5")

plt.legend()

plt.show()

