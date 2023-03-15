import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 250
nx = 64
ny = 64
step_arr = np.arange(0, ns*51, ns)
for step in step_arr:
    dfc = pd.read_csv(f"data/coni/2d{step}.csv", header=None)
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

# nx = 32
# ny = 64
# step = 10000
# x= np.linspace(0, 20*ny, num=ny)
# dfc = pd.read_csv(f"32-5dx-fo/coni/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower', cmap ="binary", vmin=0.000, vmax=0.02)
# # plt.colorbar()
# plt.plot(x,matc[:, 0], label="W=100$d_0$", linewidth=0.5, color="black")
# # plt.show()

# nx = 64
# ny = 128
# step = 40000
# x= np.linspace(0, 10*ny, num=ny)
# dfc = pd.read_csv(f"64-5dx-fo/coni/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower')
# plt.plot(x, matc[:, 0], label="W=50$d_0$", linewidth=0.9, color="black")

# nx = 128
# ny = 256
# step = 160000
# x= np.linspace(0, 5*ny, num=ny)
# dfc = pd.read_csv(f"128-5dx-fo/coni/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower', cmap ="binary", vmin=0.000, vmax=0.02)
# # plt.colorbar()
# plt.plot(x, matc[:, 0], label="W=25$d_0$", linewidth=1.4, color="black")


# nx = 256
# ny = 512
# step = 640000
# x= np.linspace(0, 2.5*ny, num=ny)
# dfc = pd.read_csv(f"256-5dx-fo/coni/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower', cmap ="binary", vmin=0.000, vmax=0.02)
# # plt.colorbar()
# plt.plot(x, matc[:, 0], label="W=12.5$d_0$", linewidth=1.9, color="black")

# nx = 256
# ny = 512
# step = 640000
# x= np.linspace(0, 2.5*ny, num=ny)
# dfc = pd.read_csv(f"256-5dx-fo/con/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower', cmap ="binary", vmin=0.000, vmax=0.02)
# # plt.colorbar()
# plt.plot(x, matc[:, 0], label="W=12.5$d_0$", linewidth=1.9, color="black")

# plt.legend()
# plt.ylim([0.0170, 0.0176])
# plt.show()
