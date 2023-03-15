import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ke = 0.3
vd = 0.68
vdb = 2.6

x = np.linspace(0, 2.6, 100)
y = (ke*(1.0-(x/vdb)**2)+x/vd)/(1.0-(x/vdb)**2+x/vd)

plt.plot(x, y, color="k", linestyle="-", linewidth=0.5, label="")

# plt.plot([0.88, 1.38, 1.78], [0.5998, 0.7402, 0.8585],
#          'o', mfc='none', markersize=8, color="k", label="$\eta$=6$\Delta$x")
# plt.plot([0.900, 1.45, 1.86], [0.5998, 0.7398, 0.8629],
#          's', mfc='none', markersize=8, color="k", label="$\eta$=12$\Delta$x")

plt.xlabel("Interface Velocity v (m/s)", fontsize=18)
plt.ylabel("Partition Coefficient", fontsize=16)
# plt.legend(fontsize=12)

# Display the plot
plt.show()


# print(dfc.loc[1, 2])
