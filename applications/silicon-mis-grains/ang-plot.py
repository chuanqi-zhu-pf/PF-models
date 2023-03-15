import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

nx = 200
ny = 150

dfp = pd.read_csv(f"data/phi/2d500.csv", header=None)
arrp = dfp[0].values
matp = arrp.reshape(nx, ny)

dfa = pd.read_csv(f"data/mob/2d500.csv", header=None)
arra = dfa[0].values
mata = arra.reshape(nx, ny)

tm = 54.7

img = np. zeros((nx, ny, 3), dtype=np. uint8)
for i in range(nx):
    for j in range(ny):
        if(matp[i, j] * (1-matp[i, j]) == 0):
            img[i, j, 0] = 255
            img[i, j, 1] = 255
            img[i, j, 2] = 255

        else:
            img[i, j, 0] = 255 * mata[i, j]/tm*(1-matp[i, j] * (1-matp[i, j]))
            img[i, j, 1] = 0
            img[i, j, 2] = 255 * (tm-mata[i, j])/tm * \
                (1-matp[i, j] * (1-matp[i, j]))

plt.imshow(img)
plt.show()
