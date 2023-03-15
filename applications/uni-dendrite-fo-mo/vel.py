import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"50d0-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]*1.0e3

plt.plot(x, y, label="$w_0$", color="black", linewidth=1)

df = pd.read_csv(f"25d0-5dx/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]*1.0e3

plt.plot(x, y, linestyle=((0, (5, 10))),
         label="$w_1$", color="black", linewidth=1.2)

# df = pd.read_csv(f"256-256-w2/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y, linestyle=(0,(1,1)),label = "$w_2$", color="black", linewidth=1)

# df = pd.read_csv(f"256-5dx-fo/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y, linestyle=(5,(10,3)), label = "w/$d_0$=12.5", linewidth=1, color="green")

# plt.ylim([0, 1.0])
plt.legend(fontsize=16)

plt.show()
