import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

N = 100
X, Y = np.mgrid[0:3:complex(0, N), 0:2:complex(0, N)]
Z1 = (1 + np.sin(Y * 10.)) * X**2

pcm = plt.pcolormesh(X, Y, Z1, norm=colors.PowerNorm(gamma=0.5))
plt.colorbar(pcm)

plt.show()
