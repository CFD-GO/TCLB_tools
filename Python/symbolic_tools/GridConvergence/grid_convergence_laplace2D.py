from sympy.abc import x
from Benchmarks.ADE.Laplace_2D_analytical import analytical_laplace_2d, InputForLaplace2DAnalytical, make_anal_plot

from DataIO.VTIFile import VTIFile
import os
import pwd
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

filename = 'laplace_benchmark_d2q9_VTK_P00_00100000.vti'
home = pwd.getpwuid(os.getuid()).pw_dir
filepath = os.path.join(home, 'DATA_FOR_PLOTS', 'Laplace_convergence', '64x64', filename)
vti_reader = VTIFile(filepath)

T_num = vti_reader.get("T")
U = vti_reader.get("U", vector=True)

T_num = np.delete(T_num, 0, axis=0)  # obligatory delete first row - wall bc (stops periodicity)

n_rows, n_columns = T_num.shape
T_num = np.delete(T_num, (n_rows - 1), axis=0)  # delete last row - extra heater bc?!

# -------- analytical solution ---------------
ySIZE, xSIZE = T_num.shape
step = 1
my_fun = -4 * x * (x - xSIZE) / (xSIZE * xSIZE)
anal_input = InputForLaplace2DAnalytical(xSIZE, ySIZE, step, my_fun)

dump_fname = f'T_anal_x{xSIZE}y{ySIZE}.npy'

if os.path.isfile(dump_fname):
    print(f'{dump_fname} found, loading results from disc')
    T_anal = np.load(dump_fname)
    x_grid = np.linspace(0, xSIZE, xSIZE)
    y_grid = np.linspace(0, ySIZE, ySIZE)
    xx, yy = np.meshgrid(x_grid, y_grid)
else:
    print(f'{dump_fname} not found, starting calculations')
    xx, yy, T_anal = analytical_laplace_2d(anal_input)
    np.save(dump_fname, T_anal)

# make_anal_plot(xx, yy, T_anal)
# -------- error norm ---------------
# T_err_field = (T_anal - T_num) / T_anal
# T_err_field[np.isnan(T_err_field)]=1
# T_err_field[np.isinf(T_err_field)]=1
# T_err_field = np.clip(T_err_field, -1, 1)

T_err_field = T_anal - T_num
T_norm = np.sum(np.sqrt((T_anal - T_num) * (T_anal - T_num)))  # Euclidean norm
print(f"T Euclidean norm={T_norm}")

print("---------- PLOTTING -------------")

fig = plt.figure(figsize=(12, 8))
ax = fig.gca(projection='3d')

# alpha=1, rstride=1, cstride=1)
ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5, label='T_err_field')
# ax.plot_surface(xx, yy, T_num,  cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, label='T_num')
# ax.plot_surface(xx, yy, T_anal,  cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.1, label='T_anal')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Customize the z axis.
# ax.set_zlim(-.1, 1.05)
# ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

plt.title(f'Laplace benchmark\n '
          f'x={xSIZE}[lu] y={ySIZE}[lu] '
          r'$T_{err}$' + f'={T_norm:.2f}')
plt.grid(True)

fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
plt.show()

print("bye")
