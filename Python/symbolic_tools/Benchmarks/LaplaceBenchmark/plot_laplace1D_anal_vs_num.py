from Benchmarks.LaplaceBenchmark.Laplace_2D_analytical import prepare_anal_data_new

from DataIO.VTIFile import VTIFile
import os
import pwd

import matplotlib.pyplot as plt
import numpy as np

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

lattice_size = 32
filename_vtk = f'laplace_template_nx_{lattice_size}_ny_{lattice_size + 2}_VTK_P00_00250000.vti'

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'LaplaceBenchmark')
# folder = os.path.join(main_folder, 'eq_sin_scheme_laplace_template')
folder = os.path.join(main_folder, 'abb_sin_scheme_laplace_template')


filepath_vtk = os.path.join(folder, filename_vtk)

vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")
U = vti_reader.get("U", is_vector=True)

# ---------------------- clip buffer bc  --------------------
T_num = np.delete(T_num, 0, axis=0)  # obligatory delete first row - wall bc (stops periodicity)

n_rows, n_columns = T_num.shape
T_num = np.delete(T_num, (n_rows - 1), axis=0)  # delete last row - extra heater bc

# ---------------------- calculate solution --------------------

xx, yy, T_anal = prepare_anal_data_new(*T_num.shape, folder, shall_recalculate_results=True)
# ---------------------- clip again --------------------
# T_num = peel_the_skin(T_num)
# T_anal = peel_the_skin(T_anal)
# xx = peel_the_skin(xx)
# yy = peel_the_skin(yy)

ySIZE, xSIZE = T_num.shape


T_mse = np.sum((T_anal - T_num) * (T_anal - T_num)) / len(T_anal)
T_L2 = np.sqrt(
    np.sum((T_anal - T_num) * (T_anal - T_num))
    / np.sum(T_anal * T_anal))  # Eq. 4.57

fig_name = f'LaplaceBenchmark_LB_vs anal_{lattice_size}x{lattice_size}.png'

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

plt.subplot(2, 1, 1)
plt.title(f'Laplace benchmark\n '
          f'lattice size={T_num.shape}[lu] '
          r'$T_{L2}$' + f'={T_L2:.2e}'
          )
plt.xlabel(r'x')
plt.ylabel(r'$Temperature$')
plt.plot(xx[xSIZE-1, :], T_anal[int(xSIZE-1), :],
         color="black", marker=">", markevery=1, markersize=3, linestyle="-", linewidth=2,
         label='analytical solution')

plt.plot(xx[xSIZE-1, :],  T_num[int(xSIZE-1), :],
         color="black", marker="o", markevery=1, markersize=3, linestyle=":", linewidth=2,
         label='current model')
axes = plt.gca()
axes.set_ylim([-0.05, 1.05])
plt.legend()
plt.grid()


plt.subplot(2, 1, 2)
plt.xlabel(r'y')
plt.ylabel(r'$Temperature$')
plt.plot(yy[:, ySIZE-1], T_anal[:, int(ySIZE/2)],
         color="black", marker=">", markevery=1, markersize=3, linestyle="-", linewidth=2,
         label='analytical solution')

plt.plot(yy[:, ySIZE-1], T_num[:, int(ySIZE/2)],
         color="black", marker="o", markevery=1, markersize=3, linestyle=":", linewidth=2,
         label='current model')
axes = plt.gca()
axes.set_ylim([-0.05, 1.05])
plt.legend()
plt.grid()


plt.show()

###################################################################################################################
# 2d clip
# xx = xx[xSIZE-1, :]
# yy = yy[:, xSIZE-1]
#
# # T_anal = T_anal[:, int(xSIZE/2)]
# T_anal = T_anal[int(xSIZE-1), :]  # dla y
#
# # T_num = T_num[:, int(xSIZE/2)]
# T_num = T_num[int(xSIZE-1), :]  # dla y
# fig_name = f'LaplaceBenchmark_LB_vs anal_{lattice_size}x{lattice_size}.png'
#
#
# plt.rcParams.update({'font.size': 14})
# plt.figure(figsize=(14, 8))
#
# axes = plt.gca()
# plt.plot(xx, T_anal,
#          color="black", marker=">", markevery=1, markersize=3, linestyle="-", linewidth=2,
#          label='analytical solution')
#
# plt.plot(yy, T_num,
#          color="black", marker="o", markevery=1, markersize=3, linestyle=":", linewidth=2,
#          label='current model')
#
#
# # ------ format y axis ------ #
# yll = T_anal.min()
# yhl = T_anal.max()
# axes.set_ylim([0, 1.0])
# # axes.set_yticks(np.linspace(yll, yhl, 5))
# # axes.set_yticks(np.arange(yll, yhl, 1E-2))
# # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# # axes.yaxis.set_major_formatter(xfmt)
#
# # plt.yscale('log')
#
# # ------ format x axis ------ #
# # plt.xlim(0, int(xSIZE / 2))
#
# # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)
#
#
# plt.title(f'Laplace benchmark\n '
#           f'lattice size={T_num.shape}[lu] '
#           r'$T_{L2}$' + f'={T_L2:.2e}'
#           )
# plt.xlabel(r'x or y')
# plt.ylabel(r'$Temperature$')
# plt.legend()
# plt.grid()
#
# fig = plt.gcf()  # get current figure
# fig.savefig(fig_name, bbox_inches='tight')
# plt.show()

# plt.close(fig)  # close the figure
