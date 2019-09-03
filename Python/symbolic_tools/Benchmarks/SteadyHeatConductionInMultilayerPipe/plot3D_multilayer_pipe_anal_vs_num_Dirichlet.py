from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeDirichlet
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import os
import pwd
from DataIO.VTIFile import VTIFile

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

lattice_size = 32
gauge = 8
# bc_scheme = 'abb_scheme'
bc_scheme = 'eq_scheme'
home = pwd.getpwuid(os.getuid()).pw_dir

# filename_vtk = f'ruraWrurze_cht_VTK_P00_00500000.vti'
# main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'ruraWrurzeBenchmark_const_k', bc_scheme, f'ruraWrurze_cht_{int(gauge * lattice_size)}lu')

k_outer = 0.01  # 0.001
filename_vtk = f'k_outer_{k_outer}_size_{int(gauge * lattice_size)}lu_VTK_P00_00500000.vti'
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_ruraWrurze_variable_k', bc_scheme, f'k_outer_{k_outer}_size_{int(gauge * lattice_size)}lu')


filepath_vtk = os.path.join(main_folder, filename_vtk)
vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")


# ----------------------- calc dimensions -----------------------

ySIZE, xSIZE = T_num.shape
assert ySIZE == xSIZE == int(gauge * lattice_size)

r0 = gauge * (8 / 2)  # inner radius
r2 = gauge * (30 / 2)  # outer radius

# abb_correction = 0.5
# if bc_scheme =='abb_scheme':
#     r0 += abb_correction
#     r2 -= abb_correction

r1 = gauge * (20 / 2)  # interface between layers

k1 = 0.1  # inner layer - heat conductivity for r0 < r < r1
k2 = 0.01  # outer layer - heat conductivity for r1 < r < r2

T0 = 0  # temperature for r = r0
T2 = 1  # temperature for r = r2

x0 = gauge * (lattice_size / 2)  # center of the pipe
y0 = gauge * (lattice_size / 2)

# ----------------------- compute anal solution ---------------------------

pwp = PipeWithinPipeDirichlet(r0, r1, r2, k1, k2, T0, T2)

x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)
T_anal = np.zeros((ySIZE, xSIZE))

for i in range(ySIZE):
    # print(f"=== Doing i/ny: {i}/{ny}  ===")
    for j in range(xSIZE):
        # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
        r = pwp.get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
        T_anal[i][j] = pwp.get_temperature_r(r)

T_err_field = T_anal - T_num
# T_L2 = np.sum(abs(T_err_field)

# -------- error norm hacks ---------------
# T_err_field = (T_anal - T_num) / T_anal
# T_err_field[np.isnan(T_err_field)]=0
# T_err_field[np.isinf(T_err_field)]=0
# T_err_field = np.clip(T_err_field, -1, 1)
# np.isnat()

T_L2 = np.sqrt(
    np.sum((T_anal - T_num) * (T_anal - T_num))
    / np.sum(T_anal * T_anal))  # Eq. 4.57

print("---------- PLOTTING -------------")
fig_name = f'pipe_within_pipe_k1{k1}_k2{k2}.png'

fig = plt.figure(figsize=(12, 8))
ax = fig.gca(projection='3d')

# alpha=1, rstride=1, cstride=1)
ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5, label='T_err_field')  # coolwarm
# ax.plot_surface(xx, yy, T_anal, cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.01, label='T_anal')
# ax.plot_surface(xx, yy, T_num, cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, alpha=1, label='T_num')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Customize the z axis.
# ax.set_zlim(-.1, .1)
# ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

plt.title(f'Pipe within pipe benchmark\n '
          f'x={xSIZE}[lu] y={ySIZE}[lu] '
          r'$T_{L2}$' + f'={T_L2:.4f}'
          )
plt.grid(True)

fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
plt.show()

print("bye")

