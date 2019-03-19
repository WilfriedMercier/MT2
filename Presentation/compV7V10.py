import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in', labelsize=14)

V4 =  (22*60+11)*1.0
V7 =  V4/(np.array([15, 14, 15, 21])*60 + np.array([0.055, 0.854, 0.432, 0.746]))
V10 =  V4/(np.array([14, 14, 13, 18])*60 + np.array([0.857, 0.857, 0.652, 0.87]))

plt.plot([1, 2, 4, 8], V10, ".", label="V10")
plt.plot([1, 2, 4, 8], V7, "x", label="V7")


plt.xlabel("Nombre de threads", size=20)
plt.ylabel("SpeedUp (par rapport à V4)", size=20)
plt.grid()
plt.legend(loc=2, fontsize=12)
plt.show()
