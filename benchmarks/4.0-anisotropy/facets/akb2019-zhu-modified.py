from email.policy import default
import numpy as np
import matplotlib.pyplot as plt

ang = np.arange(0, 54.7/180*np.pi, 0.001)
angd = ang/np.pi*180
al0 = 54.7/180*np.pi
ep = 0.02
del0 = 0.36
gamma0 = 0.44
coef = 1.0/1.62454431
aa = np.zeros(len(ang))
aap = np.zeros(len(ang))
aapp = np.zeros(len(ang))
ga = np.zeros(len(ang))
gas = np.zeros(len(ang))
# a0 = 1 + np.sqrt(np.cos(np.pi/4)**2+ep*ep) + np.sqrt(np.sin(np.pi/4)**2+ep*ep)
for idx, val in enumerate(ang):
    COS = np.sqrt(np.cos(val)**2+ep*ep)
    SIN = np.sqrt(np.sin(val)**2+ep*ep)
    aa[idx] = (1 + del0*(COS + np.tan(al0)*SIN))*coef
    aap[idx] = del0*coef*(-0.5*np.sin(2*val)/COS +
                          np.tan(al0)*np.sin(2*val)/SIN/2)
    aapp[idx] = 0.5*del0*coef*(2*np.cos(2*val)*(1/SIN-1/COS)-np.sin(2*val)
                               * (1/SIN/SIN*np.sin(2*val)/SIN/2 + 1/COS/COS*np.sin(2*val)/COS/2))
    ga[idx] = aa[idx]*aa[idx]*gamma0
    gas[idx] = (aa[idx]*aa[idx] + 2 *
                (aa[idx]*aapp[idx]+aap[idx]*aap[idx]))*gamma0

fig, ax1 = plt.subplots(figsize=(6, 4))

ax2 = ax1.twinx()
ax1.plot(angd, ga, "-", label="$\gamma(\\theta)$")
ax2.plot(angd, gas,  "--", label="$\gamma(\\theta)+\gamma{''}(\\theta)$")

# ax2.yticks(fontsize=16)
print(aa)

ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_xlabel('Angle of deviation from the <111>', fontsize=14)
ax1.set_xticklabels(['0', '$0^\circ$', '$10^\circ$',
                    '$20^\circ$', '$30^\circ$', '$40^\circ$', '$50^\circ$'])
ax1.set_ylabel('$\gamma(\\theta)$', fontsize=14)
ax2.set_ylabel("$\gamma(\\theta)+\gamma''(\\theta)$", fontsize=14)

fig.legend(bbox_to_anchor=(0.9, 0.7), loc='upper right',
           borderaxespad=4, fontsize=14)

# plt.figure(figsize=(10, 6))
plt.show()
