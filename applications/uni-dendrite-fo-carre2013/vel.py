import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"data/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]*1.0e3

plt.plot(x, y,  linestyle=(0, (1, 10)),
         label="w/$d_0$=100.0", color="black", linewidth=1)

# df = pd.read_csv(f"128-5dx/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:,0]
# y = df.loc[1:,1]/df.loc[1:,2]*1.0e3

# plt.plot(x, y,  linestyle=(0,(1,10)), label = "w/$d_0$=100.0", color="black", linewidth=1)

df = pd.read_csv(f"0_8_eta/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]*1.0e3

plt.plot(x, y, linestyle="dashed", label="$\eta$=0.8",
         color="black", linewidth=1)

df = pd.read_csv(f"1_6_eta/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]*1.0e3

plt.plot(x, y, linestyle=(0, (1, 1)),
         label="$\eta$=1.6", color="black", linewidth=1)

df = pd.read_csv(f"0_4_eta/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]*1.0e3

plt.plot(x, y, linestyle=(5, (10, 3)),
         label="$\eta$=0.4", linewidth=1, color="green")

# plt.ylim([0.2, 1.5])
plt.legend()

plt.show()
