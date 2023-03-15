import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 1000
nx = 128
ny = 256
step_arr = np.arange(0, ns*301, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    # plt.imshow(np.rot90(matc), origin='lower')
    # plt.axis('off')
    plt.clim(0.05, 0.25)
    plt.colorbar()
    # plt.plot(matc[:, 386])
    plt.savefig(f"fig/con/2d{step}")
    plt.close()

# # dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
# # arrc = dfc[0].values
# # matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower')
# # plt.colorbar()
# # # plt.plot(matc[60, :])
# # plt.savefig(f"fig/phi/2d{step}")
# # plt.close()

# # dfc = pd.read_csv(f"data/conl/2d{step}.csv", header=None)
# # arrc = dfc[0].values
# # matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower')
# # plt.colorbar()
# # # plt.plot(matc[:, 0])
# # plt.savefig(f"fig/conl/2d{step}")
# # plt.close()

    dfc = pd.read_csv(f"data/temp/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    plt.colorbar()
    # plt.plot(matc[:, 0])
    plt.savefig(f"fig/temp/2d{step}")
    plt.close()

# step = 58000
# dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.imshow(matc, origin='lower')
# # plt.imshow(np.rot90(matc), origin='lower')
# # plt.axis('off')
# plt.clim(0.07, 0.25)
# # plt.colorbar()
# # plt.plot(matc[:, 256])
# # plt.savefig(f"fig/con/2d{step}")
# plt.show()
