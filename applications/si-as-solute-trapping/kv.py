import numpy as np
import matplotlib.pyplot as plt

ke = 0.5
me = -100.0
vd = 0.68/10
vdb = 2.6/10

varr = np.arange(0, vdb, 0.00001)
karr = np.zeros(len(varr))
marr = np.zeros(len(varr))

for i, v in enumerate(varr):
    karr[i] = (ke*(1.0-(v/vdb)**2)+v/vd)/(1.0-(v/vdb)**2+v/vd)
    marr[i] = me*(1.0+(ke-karr[i]+karr[i]*np.log(karr[i]/ke))/(1.0-ke))
    if(np.abs(karr[i] - 0.833) < 0.0001):
        print(marr[i])


plt.plot(varr, karr)
plt.xlabel("Interface Velocity (v) / $m/s$")
plt.ylabel("Partition Coefficient (k)")
plt.show()
