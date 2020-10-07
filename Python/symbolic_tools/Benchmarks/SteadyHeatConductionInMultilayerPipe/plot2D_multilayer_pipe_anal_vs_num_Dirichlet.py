from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeDirichlet
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

reference_lattice_size = 32
gauge = 8
lattice_size = int(gauge * reference_lattice_size)


home = pwd.getpwuid(os.getuid()).pw_dir

k_inner = 0.1  # inner layer - heat conductivity for r0 < r < r1
k_outer = 0.1  # outer layer - heat conductivity for r1 < r < r2

# bc_scheme = 'abb_scheme'
bc_scheme = 'eq_scheme'

filename_vtk = f'k_outer_{k_outer}_size_{int(gauge * reference_lattice_size)}lu_VTK_P00_00500000.vti'
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_ruraWrurze_variable_k', bc_scheme, f'k_outer_{k_outer}_size_{lattice_size}lu')
filepath_vtk = os.path.join(main_folder, filename_vtk)
vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")


# ----------------------- calc dimensions -----------------------

ySIZE, xSIZE = T_num.shape
assert ySIZE == xSIZE == int(gauge * reference_lattice_size)

r0 = gauge * (8 / 2)  # inner radius
r2 = gauge * (30 / 2)  # outer radius

abb_correction = 0.5
if bc_scheme =='abb_scheme':
    r0 += abb_correction
    r2 -= abb_correction

r1 = gauge * (20 / 2)  # interface between layers

x0 = gauge * (reference_lattice_size / 2)  # center of the pipe
y0 = gauge * (reference_lattice_size / 2)


# ----------------------- compute anal solution ---------------------------

pwp = PipeWithinPipeDirichlet(r0, r1, r2, k_inner, k_outer, T0=0, T2=1)

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


T_err_field_eq = T_anal - T_num
# T_L2 = np.sum(abs(T_err_field)

# -------- error norm hacks ---------------
# T_err_field = (T_anal - T_num) / T_anal
# T_err_field[np.isnan(T_err_field)]=0
# T_err_field[np.isinf(T_err_field)]=0
# T_err_field = np.clip(T_err_field, -1, 1)
# np.isnat()

# 2d clip
T_anal = T_anal[:, int(xSIZE / 2)]
T_num = T_num[:, int(xSIZE / 2)]
T_mse_eq = np.sum((T_anal - T_num) * (T_anal - T_num)) / len(T_anal)
T_L2_eq = np.sqrt(
    np.sum((T_anal - T_num) * (T_anal - T_num))
    / np.sum(T_anal * T_anal))  # Eq. 4.57

###################################################################################################################

fig_name = f'pipe_within_pipe_k1_{k_inner}_k2_{k_outer}_{xSIZE}lu.png'


# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()
plt.plot(x_grid, T_anal,
         color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='analytical solution')

plt.plot(x_grid, T_num,
         color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
         label='current model')


# ------ format y axis ------ #
yll = T_anal.min()
yhl = T_anal.max()
axes.set_ylim([yll, 1.05*yhl])
# axes.set_yticks(np.linspace(yll, yhl, 5))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
# ------ format x axis ------ #
plt.xlim(0, int(xSIZE / 2))

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'T across pipes, {bc_scheme} \n '
          r'$k_{inner}$=' + f'{k_inner} \t;\t' 
          r'$k_{outer}$=' + f'{k_outer} \n'
          f'{xSIZE}x{xSIZE} [lu]'
          f'\t' r'$T_{MSE}$=' + f'{T_mse_eq:.2e}'
          )

plt.xlabel(r'$r$')
plt.ylabel(r'$Temperature$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure
