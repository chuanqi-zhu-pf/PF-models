import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.figure(figsize=(12, 9), dpi=80)
dfc = pd.read_csv(f"1590-03092/con/1d60000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black", label="1593 K", linewidth=1.5)
dfc = pd.read_csv(f"1590-03092/conl/1d60000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black",
         linewidth=1.5, linestyle="dashed")
dfc = pd.read_csv(f"1580-03092/con/1d60000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black", label="1583 K", linewidth=1.0)
dfc = pd.read_csv(f"1580-03092/conl/1d60000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black",
         linewidth=1.0, linestyle="dashed")
dfc = pd.read_csv(f"1570-03092/con/1d60000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black", label="1573 K",  linewidth=0.5)
dfc = pd.read_csv(f"1570-03092/conl/1d60000.csv", header=None)
arrc = dfc[0].values
plt.plot(arrc[100:200], color="black",
         linewidth=0.5, linestyle="dashed")
plt.ylim([0.087, 0.160])
plt.xlabel("Position ($\\times 0.5 $nm)", fontsize=14)
plt.ylabel("As Concentration (at.)", fontsize=14)
plt.legend(fontsize=14)
plt.show()
