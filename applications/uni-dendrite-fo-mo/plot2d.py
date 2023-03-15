import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 4000
nx = 128
ny = 128
step_arr = np.arange(0, ns*21, ns)
for step in step_arr:
    dfc = pd.read_csv(f"25d0-5dx/coni/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    # plt.imshow(matc, origin='lower')
    plt.plot(matc[:, 0], marker="o")
    dfc = pd.read_csv(f"25d0-5dx/con/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    # plt.imshow(matc, origin='lower')
    plt.plot(matc[:, 0], marker="o")
    plt.savefig(f"figures/con/2d{step}")
    plt.close()

    dfc = pd.read_csv(f"25d0-5dx/phi/2d{step}.csv", header=None)
    arrc = dfc[0].values
    matc = arrc.reshape(ny, nx)
    plt.imshow(matc, origin='lower', cmap='binary')
    # plt.axis('off')
    # plt.plot(matc[:, 0])
    plt.savefig(f"figures/phi/2d{step}")
    plt.close()

# nx = 64
# ny = 64
# step = 5000
# x = np.linspace(0, 20*ny, num=ny)
# dfc = pd.read_csv(f"64-5dx/con/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.plot(x, matc[:, 0],  marker="s",markersize=4, color="black", linewidth=0.5, label="$c_L$, $W_1$")
# dfc = pd.read_csv(f"64-5dx/coni/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower', cmap="binary", vmin=0.000, vmax=0.02)
# # plt.colorbar()
# plt.plot(x, matc[:, 0], marker="s",markersize=4, color="black", linewidth=1.0, label="$c_i$, $W_1$")
# # plt.show()

# nx = 128
# ny = 128
# step = 20000
# x = np.linspace(0, 10*ny, num=ny)
# dfc = pd.read_csv(f"128-5dx/con/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.plot(x, matc[:, 0], marker="o",markersize=4, linestyle="dashed", color="black", linewidth=0.5, label="$c_L$, $W_2$")
# dfc = pd.read_csv(f"128-5dx/coni/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# # plt.imshow(matc, origin='lower', cmap="binary", vmin=0.000, vmax=0.02)
# # plt.colorbar()
# # plt.hlines(matc[:, 0], xmin, xmax, colors=None, linestyles='solid')
# plt.plot(x, matc[:, 0], marker="o",markersize=4, linestyle="dashed", color="black", linewidth=1.0, label="$c_i$, $W_2$")
# plt.xlim([960, 1080])
# plt.ylim([0.0168, 0.0182])
# # plt.legend(fontsize=16)
# plt.legend(bbox_to_anchor=(0.86, 1), loc='upper left', borderaxespad=0, fontsize=15)
# plt.show()

# nx = 256
# ny = 256
# step = 40000
# x= np.linspace(0, 5*ny, num=ny)
# dfc = pd.read_csv(f"256-5dx/phi/2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.imshow(matc, origin='lower', cmap ="binary")
# plt.colorbar()
# # plt.plot(x, matc[:, 0])

# nx = 512
# ny = 512
# step = 160000
# # x= np.linspace(0, 5*ny, num=ny)
# dfc = pd.read_csv(f"2d{step}.csv", header=None)
# arrc = dfc[0].values
# matc = arrc.reshape(ny, nx)
# plt.imshow(matc, origin='lower', cmap ="binary")#, vmin=0.000, vmax=0.02)
# # plt.colorbar()
# # plt.plot(x, matc[:, 0])

# plt.show()
