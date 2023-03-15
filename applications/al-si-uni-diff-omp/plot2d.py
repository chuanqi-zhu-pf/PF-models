import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 20000
nx = 128
ny = 256
step_arr = np.arange(0, ns*21, ns)
# step_arr = [8*ns, 9*ns, 10*ns, 11*ns, 12*ns]
for step in step_arr:
    dfc = pd.read_csv(f"data/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower')
    # plt.plot(matc[:, 0])
    plt.savefig(f"figures/con/2d{step}")
    plt.close()

    dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower', cmap='binary')
    # plt.axis('off')
    # plt.plot(matc[:, 0])
    plt.savefig(f"figures/phi/2d{step}")
    plt.close()

# nx = 512
# ny = 512
# step = 110000
# dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
# arrc = dfc[0].values
# # plt.figure(figsize=(8, 8), dpi=64)
# # xx = []
# # yy = []
# matc = arrc.reshape(ny, nx)
# # for i in range(nx):
# #     for j in range(ny):
# #         if(matc[i, j] == 0.0):
# #             xx.append(i)
# #             yy.append(j)

# # plt.scatter(xx, yy, color="grey", marker=".")
# plt.imshow(matc, cmap="gray", origin="lower")
# # plt.plot(matc[:, 64])
# plt.show()
