import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

fig = plt.figure()
ax = Axes3D(fig)

Z = np.load("DEM.npy")

X, Y = np.mgrid[:Z.shape[0],:Z.shape[1]]
ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)

ax.set_xlabel('X dirn')
ax.set_ylabel('Y dirn')
ax.set_zlabel('Elevation')

plt.show()
