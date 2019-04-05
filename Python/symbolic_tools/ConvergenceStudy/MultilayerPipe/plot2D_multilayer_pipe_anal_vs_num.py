from Benchmarks.ADE.steady_two_layer_cylinder_analytical_2D import InputForMultiLayeredPipe, PipeWithinPipe
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pwd
from DataIO.VTIFile import VTIFile

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

lattice_size = 128
gauge = 1
filename_vtk = f'ruraWrurze_VTK_P00_00050000.vti'

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'RuraWRurzeBenchmark', f'ruraWrurze_{lattice_size}lu')
filepath_vtk = os.path.join(main_folder, filename_vtk)
# filepath_vtk = os.path.join(main_folder, 'DATA_FOR_PLOTS', 'abb_laplace_template', filename_vtk)

vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")

# import pandas as pd
# filepathU = os.path.join(main_folder, f'ruraWrurze_TXT_P00_00050000_T.txt')
# dataT = pd.read_csv(filepathU, delimiter=" ")  # TODO bad shape? - LLW

# ----------------------- calc dimensions -----------------------

ySIZE, xSIZE = T_num.shape
assert ySIZE == xSIZE == int(gauge * lattice_size)

r0 = gauge * (32 / 2)  # inner radius
r2 = gauge * (126 / 2)  # outer radius
r1 = (r0 + r2) / 2  # interface between layers

k1 = 1  # inner layer - heat conductivity for r0 < r < r1
k2 = 1  # outer layer - heat conductivity for r1 < r < r2

T0 = 0  # temperature for r = r0
T2 = 1  # temperature for r = r2

x0 = gauge * (lattice_size / 2)  # center of the pipe
y0 = gauge * (lattice_size / 2)

step = 1

# ----------------------- compute anal solution ---------------------------

anal_input = InputForMultiLayeredPipe(r0, r1, r2, k1, k2, T0, T2)
pwp = PipeWithinPipe(anal_input)

nx = int(xSIZE / step)
ny = int(ySIZE / step)
x_grid = np.linspace(0, xSIZE, nx)
y_grid = np.linspace(0, ySIZE, ny)
xx, yy = np.meshgrid(x_grid, y_grid)
T_anal = np.zeros((ny, nx))

for i in range(ny):
    # print(f"=== Doing i/ny: {i}/{ny}  ===")
    for j in range(nx):
        # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
        T_anal[i][j] = pwp.get_temperature_xy(xx[i][j], yy[i][j], x0, y0)

T_err_field = T_anal - T_num

print("---------- PLOTTING -------------")
fig_name = f'pipe_within_pipe_k1{k1}_k2{k2}.png'

fig = plt.figure(figsize=(12, 8))
ax = fig.gca(projection='3d')

# alpha=1, rstride=1, cstride=1)
ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5,
                label='T_err_field')  # coolwarm
ax.plot_surface(xx, yy, T_num, cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, label='T_num')
ax.plot_surface(xx, yy, T_anal, cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.1, label='T_anal')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Customize the z axis.
# ax.set_zlim(-.1, 1.05)
# ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

plt.title(f'Pipe within pipe benchmark\n '
          f'x={xSIZE}[lu] y={ySIZE}[lu] '
          # r'$T_{mse}$' + f'={T_mse:.4f}'
          )
plt.grid(True)

fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
plt.show()

print("bye")
