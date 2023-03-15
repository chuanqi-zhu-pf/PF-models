import matplotlib.pyplot as plt
import numpy as np

# (hkl)
v1 = np.array([1, -1, -10])
v1 = v1/np.linalg.norm(v1)
# ND
n1 = np.array([0, 0, 1])
# rotation axis for aligning with the ND
ax1 = np.cross(v1, n1)
ax1 = ax1/np.linalg.norm(ax1)
# rotation angle
th1 = np.arccos(np.dot(v1, n1))

# (hkl) => // ND
v1p = v1*np.cos(th1) + np.cross(ax1, v1)*np.sin(th1) + \
    ax1*np.dot(ax1, v1)*(1-np.cos(th1))

# [uvw]
v2 = np.array([-4, 6, -1])
v2 = v2/np.linalg.norm(v2)
# RD
n2 = np.array([1, 0, 0])

# [uvw] => ^ ND
v2p = v2*np.cos(th1) + np.cross(ax1, v2)*np.sin(th1) + \
    ax1*np.dot(ax1, v2)*(1-np.cos(th1))

# solution for [uvw] // RD
ax2 = np.cross(v2p, n2)
ax2 = ax2/np.linalg.norm(ax2)
th2 = np.arccos(np.dot(v2p, n2))

print("Rotation axis for ND: ", ax1[0], ax1[1], ax1[2], "\n")
print("Rotation angle for ND: ", th1/np.pi*180, "\n")
print("Rotation axis for RD: ", ax2[0], ax2[1], ax2[2], "\n")
print("Rotation angle for RD: ", th2/np.pi*180, "\n")

# verify the solution
V1 = np.array([1, -1, -10])
V1 = V1/np.linalg.norm(V1)
V2 = np.array([-4, 6, -1])
V2 = V2/np.linalg.norm(V2)
Ax1 = np.array([-0.7071, -0.7071, 0.0])
Ax1 = Ax1/np.linalg.norm(Ax1)
Ax2 = np.array([0.0, 0.0, 1.0])
Th1 = 171.95
Th1 = Th1/180*np.pi
Th2 = 33.799
Th2 = Th2/180*np.pi

V1p = V1*np.cos(Th1) + np.cross(Ax1, V1)*np.sin(Th1) + \
    Ax1*np.dot(Ax1, V1)*(1-np.cos(Th1))

V2p = V2*np.cos(Th1) + np.cross(Ax1, V2)*np.sin(Th1) + \
    Ax1*np.dot(Ax1, V2)*(1-np.cos(Th1))

V1pp = V1p*np.cos(Th2) + np.cross(Ax2, V1p)*np.sin(Th2) + \
    Ax2*np.dot(Ax2, V1p)*(1-np.cos(Th2))

V2pp = V2p*np.cos(Th2) + np.cross(Ax2, V2p)*np.sin(Th2) + \
    Ax2*np.dot(Ax2, V2p)*(1-np.cos(Th2))

print(V1pp)
print(V2pp)
print(np.dot(V1pp, V2pp))
