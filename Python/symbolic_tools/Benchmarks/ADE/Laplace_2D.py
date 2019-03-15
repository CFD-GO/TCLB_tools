# References:
# https://www.wire.tu-bs.de/lehre/ws15/pde1/lecture_2.pdf
# https://www.wire.tu-bs.de/lehre/ws15/pde1/lecture_3.pdf


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from sympy import exp, pi, integrate

import numpy as np
import sympy as sp
from sympy.abc import x, y, t
from sympy import pretty_print
import time
from sympy import fourier_series, pi, cos, sin

my_fun = x * (1 - x)  # -x*x +x

x_low = 0
x_high = 1
y_low = 0
y_high = 1
step = 0.05

n_fourier_terms = 25

lim = (x, x_low, x_high)
L = x_high - x_low
w0 = pi / L
u_sol = 0
c = [0]

start = time.process_time()
print("---------- Calculating Fourier coeff -------------")
for k in range(1, n_fourier_terms):  # skip zero
    result = (2 / L) * integrate(my_fun * sin(k * w0 * x), lim)
    c.append(result / (sp.sinh(k * pi)))
    u_sol += c[k] * sp.sin(k * pi * x) * sp.sinh(k * pi * y)

print("---------- Calculating Values -------------")

nx = int((x_high - x_low) / step)
ny = int((y_high - y_low) / step)
x_grid = np.linspace(x_low, x_high, nx)
y_grid = np.linspace(y_low, y_high, ny)
X, Y = np.meshgrid(x_grid, y_grid)
Z = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        Z[i][j] = u_sol.subs({'x': X[i][j], 'y': Y[i][j]})

print(f'\n\n Done in {time.process_time() - start} [s].')
print("---------- PLOTTING -------------")

fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=1, antialiased=False)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Customize the z axis.
# ax.set_zlim(-1.01, 1.01)
# ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

print("bye")
