import numpy as np
import matplotlib.pyplot as plt

ke = 0.3
Tm = 1685
me = -400
T = 1500

carr = np.arange(0, 0.12, 0.001)
fsarr = np.zeros(len(carr))
flarr = np.zeros(len(carr))

for i, c in enumerate(carr):
    fsarr[i] = c*np.log(c)+(1-c)*np.log(1-c) - c*np.log(ke) + \
        (1-c)*np.log((1+(Tm-T)/me)/(1+ke*(Tm-T)/me))
    flarr[i] = c*np.log(c)+(1-c)*np.log(1-c)


plt.plot(carr, flarr)
plt.plot(carr, fsarr)
plt.show()
