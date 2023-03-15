import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 300
nx = 300
ny = 300
step_arr = np.arange(0, ns*11, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    # plt.plot(matc[:, 24])
    plt.savefig(f"fig/2d{step}")
    plt.close()
