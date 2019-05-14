from sympy.abc import x
from Benchmarks.ADE.Laplace_2D_analytical import analytical_laplace_2d, InputForLaplace2DAnalytical
from Benchmarks.ADE.Laplace_2D_analytical import prepare_anal_data

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

lattice_size = 32
filename_vtk = f'laplace_template_nx_{lattice_size}_ny_{lattice_size + 2}_VTK_P00_00250000.vti'

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'LaplaceBenchmark')
folder = os.path.join(main_folder, 'eq_scheme_laplace_template')
# folder = os.path.join(main_folder, 'abb_laplace_template')
filepath_vtk = os.path.join(folder, filename_vtk)

vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")
U = vti_reader.get("U", vector=True)

T_num = np.delete(T_num, 0, axis=0)  # obligatory delete first row - wall bc (stops periodicity)

n_rows, n_columns = T_num.shape
T_num = np.delete(T_num, (n_rows - 1), axis=0)  # delete last row - extra heater bc?!

xx, yy, T_anal = prepare_anal_data(*T_num.shape, folder)


# make_anal_plot(xx, yy, T_anal)
# -------- error norm ---------------
# T_err_field = (T_anal - T_num) / T_anal
# T_err_field[np.isnan(T_err_field)]=1
# T_err_field[np.isinf(T_err_field)]=1
# T_err_field = np.clip(T_err_field, -1, 1)

T_err_field = T_anal - T_num
T_mse = np.sum((T_anal - T_num) * (T_anal - T_num))/len(T_anal)
print(f"T_mse={T_mse}")

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
          f'lattice size={T_num.shape}[lu] '
          # r'$T_{mse}$' + f'={T_mse:.4f}'
          )
plt.grid(True)

fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
# ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
plt.show()

print("bye")
