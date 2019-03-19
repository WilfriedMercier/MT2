import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in')


#Data
T, X, Y = np.genfromtxt("data/ifortdata", unpack=True)

uniq = np.unique(T)
color = cm.rainbow(np.linspace(0,1,np.size(uniq)))

for i, c in zip(uniq, color):
   print(i)
   plt.plot(X[T==i], Y[T==i], linestyle='--', color=c, label="t="+str(i))


plt.xlabel("Position (m)", size=20)
plt.ylabel("Vitesse (m/s)", size=20)
plt.grid()
plt.legend(loc=1)
plt.show()


