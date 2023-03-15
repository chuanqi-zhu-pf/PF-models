import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 8))

df = pd.read_csv(f"data/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]

plt.plot(x, y)

# df = pd.read_csv(f"64-2d/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:, 0]
# y = df.loc[1:, 1]/df.loc[1:, 2]

# plt.plot(x, y,
#          label="$\eta_a=3.00 \  \mu m$")

# df = pd.read_csv(f"128-2d/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:, 0]
# y = df.loc[1:, 1]/df.loc[1:, 2]

# plt.plot(x, y,
#          label="$\eta_b=1.50 \ \mu m$")

# df = pd.read_csv(f"256-2d/tip_pos.csv", delimiter="   ", header=None)
# x = df.loc[1:, 0]
# y = df.loc[1:, 1]/df.loc[1:, 2]

# plt.plot(x, y,
#          label="$\eta_c=0.75 \ \mu m$")

# plt.xlim([0.0, 0.034])
# plt.ylim([0.0005, 0.00125])
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.legend(fontsize=20)
plt.show()
