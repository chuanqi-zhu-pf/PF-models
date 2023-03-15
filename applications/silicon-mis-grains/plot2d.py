import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

# df = pd.read_csv(f"data/domain.csv", header=None)
ns = 500
nx = 100
ny = 200
step_arr = np.arange(0, ns*71, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')  # , cmap='binary')
    # plt.plot(matc[:, 24])
    plt.savefig(f"fig/phi/2d{step}")
    plt.close()

    dfm = pd.read_csv(f"data/mob/2d{step}.csv", header=None)
    arrm = dfm[0].values
    matm = arrm.reshape(ny, nx)
    plt.imshow(matm, origin='lower', norm=colors.PowerNorm(gamma=10.0))
    # plt.plot(matc[:, 24])
    plt.colorbar()
    plt.savefig(f"fig/mob/2d{step}")
    plt.close()
# ,
