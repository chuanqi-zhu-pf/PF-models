import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

step = 3000
plt.figure(figsize=(10, 8))
df = pd.read_csv(f"data/con/1d{step}.csv", header=None)
arr = df[0].values
plt.plot(arr, marker="o", markersize=6,
         label="$c$", color="black", linewidth=1.5)

df = pd.read_csv(f"data/conl/1d{step}.csv", header=None)
arr = df[0].values
plt.plot(arr, marker="o", mfc='none', markersize=6,
         label="$c_L$", color="black", linewidth=0.8)

df = pd.read_csv(f"data/cons/1d{step}.csv", header=None)
arr = df[0].values
plt.plot(arr, marker="o", mfc='none', markersize=6,  label="$c_S$",
         color="black", linestyle="dashed", linewidth=0.8)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=20)
plt.show()
