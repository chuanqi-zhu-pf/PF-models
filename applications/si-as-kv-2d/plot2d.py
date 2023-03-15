import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 1000
nx = 64
ny = 128
step_arr = np.arange(0, ns*21, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    # plt.plot(matc[:, 0])
    plt.colorbar()
    plt.savefig(f"fig/con/2d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    # plt.imshow(matc, origin='lower')
    plt.plot(matc[:, 0])
    plt.savefig(f"fig/phi/2d{step}")
    plt.close()
