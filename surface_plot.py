import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def plot(TwoDarray):
    """
    Plot 2-D array in xyz space using values and 2-D grid
    Input:
      TwoDarray: 2-Dimensional array 
    Output:
      3-D plot
    """
    fig = plt.figure()
    ax = Axes3D(fig)

    Z = TwoDarray
    X, Y = np.mgrid[:Z.shape[0],:Z.shape[1]]
    ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)

    ax.set_xlabel('X dirn')
    ax.set_ylabel('Y dirn')
    ax.set_zlabel('Elevation')
    plt.show()
