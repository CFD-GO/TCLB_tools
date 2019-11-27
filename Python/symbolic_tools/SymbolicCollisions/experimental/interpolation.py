import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy

# x = np.linspace(0, 4, 13)
# y = np.array([0, 2, 3, 3.5, 3.75, 3.875, 3.9375, 4])
# X, Y = np.meshgrid(x, y)
# Z = np.sin(np.pi*X/2) * np.exp(Y/2)
#
# x2 = np.linspace(0, 4, 65)
# y2 = np.linspace(0, 4, 65)
# f = interp2d(x, y, Z, kind='cubic')
# Z2 = f(x2, y2)
#
# fig, ax = plt.subplots(nrows=1, ncols=2)
# ax[0].pcolormesh(X, Y, Z)
#
# X2, Y2 = np.meshgrid(x2, y2)
# ax[1].pcolormesh(X2, Y2, Z2)
#
# plt.show()

print('bye')
# # Suppose we want to interpolate the 2-D function


def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
# on a grid in [0, 1]x[0, 1]

# >>>
grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
# but we only know its values at 1000 data points:

# >>>
points = np.random.rand(1000, 2)
values = func(points[:,0], points[:,1])
# This can be done with griddata – below we try out all of the interpolation methods:

# >>>

grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
# # One can see that the exact result is reproduced by all of the methods to some degree, but for this smooth function the piecewise cubic interpolant gives the best results:
#
# # >>>
#
plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()