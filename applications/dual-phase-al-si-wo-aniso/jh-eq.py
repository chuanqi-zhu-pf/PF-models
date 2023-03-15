import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dC = 0.984
fa = 0.1077
fb = 0.8922
tha = 3/8*np.pi
thb = 3/8*np.pi

m1 = 680
m2 = 1302
mm = m1*m2/(m1+m2)
Diff = 0.1
gta = 430
gtb = 430
vv = 1.6e-3
lf = 1.0
tf = 1.0

# m1 = 680
# m2 = 1302
# mm = m1*m2/(m1+m2)
# Diff = 2.5e-8
# gta = 1.67e-7
# gtb = 1.67e-7
# vv = 5.0e-5
# lf = 1.0
# tf = 1.0

# m1 = 680
# m2 = 1302
# mm = m1*m2/(m1+m2)
# Diff = 2.5e-8
# gta = 1.67e-7 * 1e6
# gtb = 1.67e-7 * 1e6
# vv = 5.0e-5 * 1e-2
# lf = 1e4
# tf = 1e2

ks = 0
for i in range(1, 10):
    ks += np.sin(i*np.pi*fa)*np.sin(i*np.pi*fa)/np.pi/np.pi/np.pi/i/i/i

Kc = mm*dC/fa/fb/Diff*ks
Kr = 2.0*mm*(gta*np.cos(tha)/m1/fa + gtb*np.cos(thb)/m2/fb)

print(Kr/gta/Diff/Kc)

lex = np.sqrt(Kr/Kc/vv)/2
lmin = lex/2
lmax = lex*4

sp = np.linspace(lmin, lmax, num=100)
dt = np.zeros(100)

fig, ax = plt.subplots()
dt = (Kc*vv*sp + Kr/sp)
ax.plot(sp/2 / lf, dt/tf)

plt.legend()
plt.show()

# Kc = mm*dC/2/np.pi/Diff
# Kr = 2*((fb*m2*gta*np.sin(tha) + fa*m1*gtb*np.sin(thb)))/fa/fb/(m1+m2)

# vv = np.linspace(1.0e-6, 100.0e-6, num=100)
# sp = np.sqrt(Kr/Kc/100e-6)
# plt.plot(vv, sp, color="black", linewidth=1.0)
# print(sp)

# sp = np.sqrt(Kr/Kc/500e-6)
# print(sp)

# sp = np.sqrt(Kr/Kc/900e-6)
# print(sp)

# dt = Kc*100e-6*9.23e-6 + Kr/9.23e-6
# print(dt)

# dt = Kc*500e-6*4.13e-6 + Kr/4.13e-6
# print(dt)

# dt = Kc*900e-6*3.08e-6 + Kr/3.08e-6
# print(dt)

# + 1.0e-3/0.25/0.03)
# label="100 $\mu m/s$")

# ax.ylim([0.40, 0.60])

# dt = (Kc*0.25e-3*sp + Kr/sp)
# plt.plot(sp/2, dt)  # label="500 $\mu m/s$")

# dt = (Kc*900e-6*sp + Kr/sp)
# plt.plot(sp, dt, label="900 $\mu m/s$")

# df = pd.read_csv("dt-sp.csv", header=None, sep=' ')
# sp = df.loc[0].values
# v1 = df.loc[1].values

# ax2 = ax.twinx()
# ax.plot(sp, v1, marker='+')
# plt.ylim([0.40, 0.60])
