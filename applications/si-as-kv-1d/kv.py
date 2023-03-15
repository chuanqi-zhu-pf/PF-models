import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ke = 0.3
vd = 0.68
vdb = 2.6

x = np.linspace(0, 2.6, 100)
y = (ke*(1.0-(x/vdb)**2)+x/vd)/(1.0-(x/vdb)**2+x/vd)

plt.plot(x, y, color="k", linestyle="-", linewidth=0.5, label="")

plt.plot([0.81, 1.63, 2.21], [0.696, 0.863, 0.955],
         'o', mfc='none',markersize=8, color="k", label="$\eta$=8$\Delta$x")
plt.plot([0.786, 1.65, 2.267], [0.692, 0.862, 0.953],
         's', mfc='none',markersize=8, color="k", label="$\eta$=16$\Delta$x")

plt.xlabel("Interface Velocity v (m/s)", fontsize=14)
plt.ylabel("Partition Coefficient", fontsize=12)
plt.legend(fontsize=12)

# Display the plot
plt.show()


# print(dfc.loc[1, 2])
