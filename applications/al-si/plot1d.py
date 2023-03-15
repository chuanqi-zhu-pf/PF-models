import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 10000
step_arr = np.arange(0, ns*13, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arrc = dfc[0].values
    plt.plot(arrc, color="black", linewidth=1)
    plt.savefig(f"fig/phi/1d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/con/1d{step}.csv", header=None)
    arrc = dfc[0].values
    plt.plot(arrc, color="black")
    plt.savefig(f"fig/con/1d{step}")
    plt.close()
