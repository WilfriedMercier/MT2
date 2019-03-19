import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in', labelsize=14)

V4    =  (22*60+11)*1.0
V6    =  V4/(np.array([22, 18, 13, 14, 15])*60.0 + np.array([0.17, 0.33, 0.126, 0.689, 0.936]))
V7    =  V4/(np.array([22, 13, 7, 8, 8])*60 + np.array([0.008, 0.438, 0.89, 0.938, 0.383]))
V8    =  V4/(np.array([43, 26, 19, 14, 11])*60 + np.array([0.399, 0.047, 0.021, 0.228, 909]))
V9    =  V4/(np.array([22, 13, 7, 9, 7])*60 + np.array([0.305, 0.217, 0.866, 0.176, 0.846]))
V9G   =  V4/(np.array([24, 16, 20, 30, 40])*60 + np.array([0.197, 0.795, 0.72, 0.144, 0.99]))
VBUFF =  V4/(np.array([22, 11, 6, 0])*60 + np.array([0.438, 542, 0.303, 0.77*V4/2.775]))


plt.plot([1, 2, 4, 6, 8], V6, "-.", label="V6")
plt.plot([1, 2, 4, 6, 8], V7, "-x", label="V7")
plt.plot([1, 2, 4, 6, 8], V8, ":^", label="V8")
plt.plot([1, 2, 4, 6, 8], V9, "--v", label="V9")
plt.plot([1, 2, 4, 6, 8], V9G, "--v", label="V9 GUIDED")
plt.plot([1, 2, 4, 8], VBUFF, "-o", label="VBUFF")


plt.xlabel("Nombre de threads", size=20)
plt.ylabel("SpeedUp (par rapport à V4)", size=20)
plt.grid()
plt.legend(loc=2, fontsize=12)
plt.show()
