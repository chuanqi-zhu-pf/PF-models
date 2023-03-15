import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 2000
step_arr = np.arange(0, ns*81, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/con/1d{step}.csv", header=None)
    arrc = dfc[0].values
    plt.plot(arrc)
    plt.savefig(f"figures/con/1d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arrc = dfc[0].values
    plt.plot(arrc)
    plt.savefig(f"figures/phi/1d{step}")
    plt.close()
