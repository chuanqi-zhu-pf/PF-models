import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# df = pd.read_csv(f"data/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y,  linestyle=(0,(1,10)), label = "w/$d_0$=100.0", color="black", linewidth=1)

df = pd.read_csv(f"128-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y,  linestyle=(0,(1,10)), label = "$w_1$", color="black", linewidth=1)

df = pd.read_csv(f"256-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y, linestyle="dashed", label = "$w_2$", color="black", linewidth=1)

df = pd.read_csv(f"512-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y, linestyle=(0,(1,1)),label = "$w_3$", color="black", linewidth=1)

df = pd.read_csv(f"1024-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y, linestyle=(5,(10,3)), label = "$w_4$", linewidth=1, color="green")

plt.ylim([0.3, 1.0])
plt.legend(fontsize=16)

plt.show()

