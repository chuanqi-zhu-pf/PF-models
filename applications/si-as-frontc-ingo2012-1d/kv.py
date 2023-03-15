import numpy as np
import matplotlib.pyplot as plt

ke = 0.3
vd = 0.68
vdb = 2.6

varr = np.arange(0, vdb, 0.00001)
karr = np.zeros(len(varr))

for i, v in enumerate(varr):
    karr[i] = (ke*(1.0-(v/vdb)**2)+v/vd)/(1.0-(v/vdb)**2+v/vd)

plt.plot(varr, karr)
plt.xlabel("Interface Velocity (v) / $m/s$")
plt.ylabel("Partition Coefficient (k)")
plt.show()
