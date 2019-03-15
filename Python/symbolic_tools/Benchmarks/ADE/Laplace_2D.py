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

from sympy import fourier_series, pi, cos, sin

# from mpmath import sinh, cosh, sin, cos

# pi = sp.pi
# myfun = x * (1 - x)  # -x*x +x
#
# x_low = 0
# x_high = 1
# y_low = 0
# y_high = 1
#
# n_fourier = 8
# s = fourier_series(myfun, (x, x_low, x_high))
# # s.sigma_approximation(n)
# ns = s.truncate(n_fourier)
#
# print(f'Fourier series of {myfun} from {x_low} to {x_high}:')
# pretty_print(s)
#
# print(f'Fourier series of {myfun} from {x_low} to {x_high} truncated to first {n_fourier} terms:')
# pretty_print(ns)
# print("---------------------")

# ns
# 4*cos(x) - cos(2*x) + 4*cos(3*x)/9 - pi**2/3
# ns.args
# (-cos(2*x), 4*cos(x), -pi**2/3, 4*cos(3*x)/9)
# ns.args[1].args[0]
# 4
# ns.args[1].args[1]
# cos(x)
# g_coeff = []
# a = []
# g2 = []

# num_pi = pi.evalf()  # numerical pi
# for k in range(0, n_fourier):
#     try:
#         g_coeff.append(ns.args[k].args[0])
#         g2.append(ns.args[k].args[1])
#     except IndexError:
#         g_coeff.append(ns.args[k])
#         g2.append(None)
#
#     a.append(g_coeff[k] / sp.sinh((k + 1) * num_pi))  # loop shall start from 1 but list indices starts from 0
#     print(f'coeff {g_coeff[k]} \t and  {g2[k]}')

# for k in range(0, n_fourier):
#     try:
#         for element in ns.args[k]:
#                 g_coeff.append(element.args[0])
#                 # g2.append(ns.args[k].args[1])
#     except TypeError:
#         g_coeff.append(ns.args[k])
#         # g2.append(None)
#
#     a.append(g_coeff[k] / sp.sinh((k + 1) * num_pi))  # loop shall start from 1 but list indices starts from 0
#     print(f'coeff {g_coeff[k]}') #' \t and  {g2[k]}')



# print(ns.subs({'x': 2}).evalf())
# u_sol = a[0]*sin(k*num_pi*x)*sinh(k*num_pi*y)
# u_sol = 0
# for k in range(0, n_fourier, 1):  # loop shall start from 1 but list indices starts from 0
#     u_sol += a[k] * sp.sin((k + 1) * num_pi * x) * sp.sinh((k + 1) * num_pi * y)


# s = fourier_series(x * (1 - x), (x, 0, 1))

myfun = x * (1 - x)  # -x*x +x

x_low = 0
x_high = 1
y_low = 0
y_high = 1

# n_fourier = 6
# s = fourier_series(myfun, (x, x_low, x_high))
# print(f'Fourier series of {myfun} from {x_low} to {x_high}:')
# pretty_print(s)
#
# ns = s.truncate(n_fourier)
# print(f'Fourier series of {myfun} from {x_low} to {x_high} truncated to first {n_fourier} terms:')
# pretty_print(ns)
#
# s_iter = s.truncate(None)  #  If n is None return an iterator.
# cosine_coeffs = []
# sine_coeffs = [0]       # there is no sine term for k = 0
# for k in range(0, n_fourier):
#     term = next(s_iter)
#     print(term)
#     cosine_coeffs.append(term)
#     if 'sin' in str(term):
#         sine_coeffs.append(term)
#
# print(sine_coeffs)
#
# sine_coeffs =[ - 1/pi - 2/3, - 1/(2*pi) - 2/3, - 1/(3*pi) - 2/3, - 1/(4*pi) - 2/3, - 1/(5*pi) - 2/3, - 1/(6*pi) - 2/3, - 1/(7*pi) - 2/3]


lim = (x, x_low, x_high)
T = x_high - x_low
w0 = pi / T

my_fun = -x * x + x

u_sol = 0

c = [0]
for k in range(1, 5):
    result = (2 / T) * integrate(my_fun * sin(k * w0 * x), lim)
    c.append(result / (sp.sinh(k * pi)))
    u_sol += c[k] * sp.sin(k * pi * x) * sp.sinh(k * pi * y)

# a = [0]
# for k in range(1, len(sine_coeffs)):
#     a.append(sine_coeffs[k] / sp.sinh(k * pi))
#     u_sol += a[k] * sp.sin(k * pi * x) * sp.sinh(k * pi * y)

print("---------- Calculating values -------------")
step = 0.05
nx = int((x_high - x_low) / step)
ny = int((y_high - y_low) / step)
x_grid = np.linspace(x_low, x_high, nx)
y_grid = np.linspace(y_low, y_high, ny)
X, Y = np.meshgrid(x_grid, y_grid)
Z = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        Z[i][j] = u_sol.subs({'x': X[i][j], 'y': Y[i][j]})

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
