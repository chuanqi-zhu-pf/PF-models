import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

# df = pd.read_csv(f"data/domain.csv", header=None)
ns = 500
nx = 200
ny = 200
step_arr = np.arange(ns, ns*71, ns)
for step in step_arr:
    dfm = pd.read_csv(f"data/mob/2d{step}.csv", header=None)
    arrm = dfm[0].values
    matm = arrm.reshape(ny, nx)
    plt.imshow(matm, norm=colors.PowerNorm(gamma=4.0))
    # plt.plot(matc[:, 24])
    plt.colorbar()
    plt.savefig(f"fig/mob/2d{step}")
    plt.close()
