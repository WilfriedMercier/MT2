import matplotlib.pyplot as plt
import numpy as np

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in')


#Data
X, Y = np.genfromtxt("../data_init", unpack=True)

plt.plot(X,Y)
plt.xlabel("Position (m)", size=20)
plt.ylabel("Vitesse (m/s)", size=20)
plt.grid()
plt.show()


