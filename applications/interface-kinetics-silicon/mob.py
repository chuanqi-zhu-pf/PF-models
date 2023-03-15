import matplotlib.pyplot as plt
import numpy as np

# dtarr = np.arange(0, 0.3, 0.01)
# v1arr = np.zeros(len(dtarr))
# v2arr = np.zeros(len(dtarr))

# for i, dt in enumerate(dtarr):
#     # v1arr[i] = 0.12*dt
#     v2arr[i] = 0.00025*(dt)**2

#     # v3arr[i] = 3.42*np.exp(-14.64/dt)*(dt)**(2/3)

# # plt.plot(dtarr, v1arr)
# plt.plot(dtarr, v2arr)
# plt.ylabel("Interface Velocity v / ($m/s$)")
# plt.xlabel("Undercooling $\Delta T$ / (K)")
# plt.show()

delta = 0.001
w = 90/10
# tharr = np.arange(0, 25, 0.1)
# muarr = np.zeros(len(tharr))

# for i, th in enumerate(tharr):
#     if(th < 10):
#         muarr[i] = delta + (1-delta)*np.tan(w*th/180*np.pi) * \
#             np.tanh(1/np.tan(w*th/180*np.pi))
#     else:
#         muarr[i] = 1.0

# plt.plot(tharr, muarr)
# plt.show()

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection='3d')

x0 = np.arange(0, 25, 0.1)
x = np.zeros(len(x0))
for i, xx in enumerate(x0):
    if(xx < 10):
        x[i] = delta + (1-delta)*np.tan(w*xx/180*np.pi) * \
            np.tanh(1/np.tan(w*xx/180*np.pi))
    else:
        x[i] = 1.0
y = np.arange(0, 0.3, 0.01)

X, Y = np.meshgrid(x, y)


# Z = X*0.12*Y + (1-X)*3.42*np.exp(-14.64/Y)*(Y)**(2/3)

Z = X*0.12 + (1-X)*0.00025*Y

X, Y = np.meshgrid(x0, y)

surf = ax.plot_surface(X, Y, Z, cmap=plt.cm.cividis)

# Set axes label
ax.set_xlabel('Angle of deviation $\Theta$', labelpad=5, fontsize=15)
ax.set_ylabel('Undercooling $\Delta T$ / K', labelpad=5, fontsize=15)
ax.set_zlabel('Attachment kinetics $\mu$ / $m/(sK)$',
              labelpad=5, fontsize=15)

# fig.colorbar(surf, shrink=0.5, aspect=8)

plt.show()
