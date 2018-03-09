import pandas as pd
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from scipy import interpolate

from scipy.interpolate import griddata

# Read 3D data points of the surface of the hull form
f = open('3d_data.txt', 'r')
data = f.readlines()

X = []
Y = []
Z = []

for i in data:
    temp = i.split(' ')
    X.append(temp[0])
    Y.append(temp[1])
    Z.append(temp[2])

x = np.asarray(X, dtype=np.float32)
y = np.asarray(Y, dtype=np.float32)
z = np.asarray(Z, dtype=np.float32)

X,Y = np.meshgrid(x,y)
Z = np.tile(z, (len(z),1))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z,  cmap='summer', rstride=1, cstride=1)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

print(z)