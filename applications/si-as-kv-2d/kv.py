import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ke = 0.3
vd = 0.68
vdb = 2.6

x = np.linspace(0, 2.6, 100)
y = (ke*(1.0-(x/vdb)**2)+x/vd)/(1.0-(x/vdb)**2+x/vd)

plt.plot(x, y, color="k", linestyle="-", linewidth=0.5)

plt.plot([0.81, 1.63, 2.21], [0.696, 0.863, 0.955],
         'o', mfc='none', color="k", label="8$\Delta$x")
plt.plot([0.786, 1.65, 2.267], [0.692, 0.862, 0.953],
         's', mfc='none', color="k", label="16$\Delta$x")

plt.xlabel("Interface Velocity v (m/s)", fontsize=12)
plt.ylabel("Partition Coefficient", fontsize=12)
plt.legend(title='Interace width')

# Display the plot
plt.show()


# print(dfc.loc[1, 2])
