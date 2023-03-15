import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 5000
nx = 351
ny = 351
step_arr = np.arange(0, ns*11, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    # plt.plot(matc[:, 64])
    plt.savefig(f"fig/phi/2d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/temp/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    # plt.plot(matc[:, 64])
    plt.savefig(f"fig/temp/2d{step}")
    plt.close()
