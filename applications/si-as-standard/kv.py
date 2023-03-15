import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ke = 0.14
vd = 2.4

x = np.linspace(0, 2.4, 100)
y = (ke+x/vd)/(1.0+x/vd)

plt.plot(x, y, color="k", linestyle="-", linewidth=0.5, label="")

plt.xlabel("Interface Velocity v (m/s)", fontsize=14)
plt.ylabel("Partition Coefficient", fontsize=12)
# plt.legend(fontsize=12)

# Display the plot
plt.show()


# print(dfc.loc[1, 2])
