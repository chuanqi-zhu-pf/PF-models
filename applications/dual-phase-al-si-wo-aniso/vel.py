import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.figure(figsize=(10, 8))
df = pd.read_csv(f"data/tip_pos.csv", delimiter="   ", header=None)
x = df.loc[1:, 0]
y = df.loc[1:, 1]/df.loc[1:, 2]

plt.plot(x, y)

# plt.xlim([0.0, 0.034])
# plt.ylim([0.3, 2.5])
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.legend(fontsize=20)
plt.show()
