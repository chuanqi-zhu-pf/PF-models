import numpy as np
import matplotlib.pyplot as plt

# plt.xlabel(r"Angle from <111> direction,  $\alpha$")
# plt.ylabel(r"Relative attachment kinetic coefficient, $\beta / \beta_{100}$")

# delta = 0.707
# ep = 0.05
# w = 90/20
# ang = np.arange(0, 45, 0.1)
# aa = np.zeros(len(ang))
# for idx, val in enumerate(ang):
#     if(np.abs(val) < 20):
#         aa[idx] = delta + (1-delta)*np.sqrt(np.tan(w*val/180*np.pi)*np.tan(w*val/180*np.pi)+ep*ep) * \
#             np.tanh(1/np.sqrt(np.tan(w*val/180*np.pi)
#                     * np.tan(w*val/180*np.pi)+ep*ep))
#     else:
#         aa[idx] = 1.0
# plt.plot(ang, aa, label="{110}")

delta = 0.35
ep = 0.05
w = 90/10
ang = np.arange(0, 54.7, 0.1)
aa = np.zeros(len(ang))
for idx, val in enumerate(ang):
    if(np.abs(val) < 10):
        aa[idx] = delta + (1-delta)*np.sqrt(np.tan(w*val/180*np.pi)*np.tan(w*val/180*np.pi)+ep*ep) * \
            np.tanh(1/np.sqrt(np.tan(w*val/180*np.pi)
                    * np.tan(w*val/180*np.pi)+ep*ep))
    else:
        aa[idx] = 1.0
# plt.figure(figsize=(6, 3))
plt.plot(ang, aa, label="{111}")
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('$b(\\theta)$', fontsize=14)
plt.xlabel('$\\theta \ (^\circ)$', fontsize=14)
# plt.legend(loc='lower right', fontsize=16)

plt.show()
