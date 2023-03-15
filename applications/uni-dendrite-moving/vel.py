import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:,0]
y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

plt.plot(x, y,  linestyle=(0,(1,10)), label = "w/$d_0$=100.0", color="black", linewidth=1)

# df = pd.read_csv(f"32-5dx-fo/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y,  linestyle=(0,(1,10)), label = "w/$d_0$=100.0", color="black", linewidth=1)

# df = pd.read_csv(f"64-5dx-fo/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y, linestyle="dashed", label = "w/$d_0$=50.0", color="black", linewidth=1)

# df = pd.read_csv(f"128-5dx-fo/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y, linestyle=(0,(1,1)),label = "w/$d_0$=25.0", color="black", linewidth=1)

# df = pd.read_csv(f"256-5dx-fo/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y, linestyle=(5,(10,3)), label = "w/$d_0$=12.5", linewidth=1, color="green")

# plt.ylim([0, 0.0055])
plt.legend()

plt.show()

