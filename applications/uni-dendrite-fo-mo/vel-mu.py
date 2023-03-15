import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"64-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y, linestyle=((0, (5, 10))), label = "$\mu$1", color="black", linewidth=1.2)

df = pd.read_csv(f"64-5dx-2/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y, linestyle=(0,(1,1)),label = "$\mu$2", color="black", linewidth=1)

df = pd.read_csv(f"64-5dx-3/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y,label = "$\mu$3", color="black", linewidth=1)

plt.ylim([0, 1.0])
plt.legend(fontsize=16)

plt.show()

