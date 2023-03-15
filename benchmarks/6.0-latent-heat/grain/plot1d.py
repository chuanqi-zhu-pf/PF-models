import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 200
step_arr = np.arange(0, ns*11, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/1d{step}.csv", header=None)
    arrc = dfc[0].values
    plt.plot(arrc)
    plt.savefig(f"fig/1d{step}")
    plt.close()

    # dfc = pd.read_csv(f"data/temp/1d{step}.csv", header=None)
    # arrc = dfc[0].values
    # plt.plot(arrc)
    # plt.savefig(f"fig/temp/2d{step}")
    # plt.close()
