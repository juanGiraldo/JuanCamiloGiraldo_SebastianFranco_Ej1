import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data.txt')
V_temp = data[:, 0]
dVx_temp = data[:, 1]
dVy_temp = data[:, 2]

V = V_temp.reshape((251, 251))
dVx = dVx_temp.reshape((251, 251))
dVy = dVy_temp.reshape((251, 251))

Y, X = np.mgrid[0:5:251j, 0:5:251j]

fig, ax = plt.subplots()
strm = ax.streamplot(X, Y, dVx, dVy, density=0.6, color='k')
plt.imshow(V, cmap='jet', extent=(0, 5, 0, 5))
plt.savefig('placas.pdf')
