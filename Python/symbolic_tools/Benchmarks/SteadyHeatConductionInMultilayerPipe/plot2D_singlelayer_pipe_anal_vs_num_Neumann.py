from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeNeumann
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pwd
from DataIO.VTIFile import VTIFile
from DataIO.helpers import find_oldest_iteration

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

reference_lattice_size = 32
gauge = 4
lattice_size = int(gauge * reference_lattice_size)

CollisionType = 'CM_HIGHER'
k = f"0.1"
home = pwd.getpwuid(os.getuid()).pw_dir
# k_0.1666666_size_128lu


main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_ruraWrurze_NeumannBC', f'{CollisionType}_k_{k}_size_{lattice_size}lu')
oldest = find_oldest_iteration(main_folder)
filename_vtk = f'{CollisionType}_k_{k}_size_{int(gauge * reference_lattice_size)}lu_VTK_P00_{oldest}.vti'

filepath_vtk = os.path.join(main_folder, filename_vtk)
vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")
T_num_slice = T_num[:, :, 1]

# ----------------------- calc dimensions -----------------------

ySIZE, xSIZE = T_num_slice.shape
assert ySIZE == xSIZE == int(gauge * reference_lattice_size)

r0 = gauge * (8 / 2)  # inner radius
r2 = gauge * (30 / 2)  # outer radius

# abb_correction = 0.5
# r0 += 0.5

# J0 = 0.22  # heat flux (dT/dr) for r = r0
J0 = 0.1  # heat flux for r = r0 --> (dT/dr) = J0/k
T2 = 0  # temperature for r = r2

# ----------------------- compute anal solution ---------------------------
x0 = gauge * (reference_lattice_size / 2)  # center of the pipe
y0 = gauge * (reference_lattice_size / 2)

pwp = PipeWithinPipeNeumann(r0, r2, J0/float(k), T2)

x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)
T_anal = np.zeros((ySIZE, xSIZE))

for i in range(ySIZE):
    for j in range(xSIZE):
        r = pwp.get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
        T_anal[i][j] = pwp.get_temperature_r(r)


T_err_field_eq = T_anal - T_num_slice
# T_L2 = np.sum(abs(T_err_field)

# -------- error norm hacks ---------------
# T_err_field = (T_anal - T_num) / T_anal
# T_err_field[np.isnan(T_err_field)]=0
# T_err_field[np.isinf(T_err_field)]=0
# T_err_field = np.clip(T_err_field, -1, 1)
# np.isnat()

# nan_mask = np.argwhere(np.isnan(T_anal))
not_nan_mask = ~np.isnan(T_anal)
T_anal_masked = T_anal[not_nan_mask]
T_num_slice_masked = T_num_slice[not_nan_mask]

T_mse_eq = np.sum((T_anal_masked - T_num_slice_masked) * (T_anal_masked - T_num_slice_masked)) / len(T_anal_masked)
T_L2_eq = np.sqrt(
    np.sum((T_anal_masked - T_num_slice_masked) * (T_anal_masked - T_num_slice_masked))
    / np.sum(T_anal_masked * T_anal_masked))  # Eq. 4.57

# T_mse_eq = np.sum((T_anal - T_num_slice) * (T_anal - T_num_slice)) / len(T_anal)
# T_L2_eq = np.sqrt(
#     np.sum((T_anal - T_num_slice) * (T_anal - T_num_slice))
#     / np.sum(T_anal * T_anal))  # Eq. 4.57

# x_slice = np.arange(0, int(xSIZE / 2), 1) + 0.5
T_num_slice = T_num_slice[:, int(xSIZE / 2)]  # take Y slice
T_num_slice = T_num_slice[int(xSIZE / 2):]  # half of it

T_anal = T_anal[:, int(xSIZE / 2)]  # take Y slice
T_anal = T_anal[int(xSIZE / 2):]  # half of it

x = x_grid[:int(xSIZE / 2)]  # half of it

# step = 0.01
# r = np.arange(r0, r2, step) + 0.5
mask = (x > r0) & (x < r2)
r = x[mask]
T_r_anal = np.array([pwp.get_temperature_r(r_) for r_ in r])
T_num_slice = T_num_slice[mask]


###################################################################################################################

fig_name = f'pipe_within_pipe_Neumann_anal_vs_num_J0{J0}_T2{T2}_{xSIZE}lu.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

plt.plot(r, T_r_anal,
         color="black", marker="v", markevery=5, markersize=7, linestyle="-", linewidth=2,
         label='analytical solution')

plt.plot(r, T_num_slice,
         color="black", marker="o", markevery=5, markersize=7, linestyle=":", linewidth=2,
         label='current model')


# ------ format y axis ------ #
yll = min(T_r_anal.min(), T_num_slice.min())
# yhl = T_r_anal.max()
axes.set_ylim([yll, 1.05])
# axes.set_yticks(np.linspace(yll, yhl, 5))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
# ------ format x axis ------ #
plt.xlim(0, int(xSIZE / 2))
# plt.xlim(int(xSIZE / 2), xSIZE)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'T across pipes,  \n '
          f'{xSIZE}x{xSIZE} [lu]'
          f'\t' r'$T_{MSE}$=' + f'{T_mse_eq:.2e}'
          f'\t' r'$T_{L2}$=' + f'{T_L2_eq:.2e}'
          )

plt.xlabel(r'$r$')
plt.ylabel(r'$Temperature$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure
