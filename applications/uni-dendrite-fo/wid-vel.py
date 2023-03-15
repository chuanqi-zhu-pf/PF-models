import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(f"tip_vel.csv", delimiter=", ", header=None)
x = 1/df.loc[:,0].values
y = df.loc[:,1].values*1.0e3
plt.plot(x, y, marker="o", color="black", linewidth=1)
plt.show()
# print(x)
# print(y)