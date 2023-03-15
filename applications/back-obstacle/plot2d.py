import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 1000
nx = 64
ny = 128
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
    plt.imshow(matc, origin='lower')#, cmap='binary')
    # plt.axis('off')
    # plt.plot(matc[:, 0])
    plt.savefig(f"figures/phi/2d{step}")
    plt.close()

# nx = 64
# ny = 128
# step = 10000
# dfc = pd.read_csv(f"data/phi/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.imshow(matc, origin='lower')
# # plt.plot(matc[:, 0])
# plt.show()

# nx = 128
# ny = 256
# step = 64000
# dfc = pd.read_csv(f"128-5dx/con/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower')
# plt.plot(matc[:, 0])
