import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 2000
nx = 64
ny = 128
step_arr = np.arange(0, ns*201, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower', cmap="cividis")
    plt.colorbar()
    # plt.plot(matc[13, :])
    plt.savefig(f"fig/con/2d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower', cmap="cividis")
    plt.colorbar()
    # plt.plot(matc[:, 64])
    plt.savefig(f"fig/phi/2d{step}")
    plt.close()

# dfc = pd.read_csv(f"data/conl/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.imshow(matc, origin='lower',  cmap="cividis")
# plt.colorbar()
# # plt.plot(matc[:, 64])
# plt.savefig(f"fig/conl/2d{step}")
# plt.close()

# dfc = pd.read_csv(f"data/con/2d{ns}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.imshow(matc, origin='lower')
# # plt.plot(matc[:, 64])
# # plt.savefig(f"fig/con/2d{500}")
# plt.show()
