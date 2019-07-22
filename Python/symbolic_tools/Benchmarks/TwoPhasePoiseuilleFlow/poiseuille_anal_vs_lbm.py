import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import os
from Benchmarks.TwoPhasePoiseuilleFlow.TwoPhasePoiseuilleAnal import TwoPhasePoiseuilleAnal, calc_gx
from Benchmarks.TwoPhasePoiseuilleFlow.TwoPhasePoiseuilleFD import TwoPhasePoiseuilleFD
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, calc_mse, calc_L2
import pwd
from DataIO.VTIFile import VTIFile


rho_h = 1
rho_l = 1
# kin_visc_h = 1/6
# kin_visc_l = 1/6

kin_visc_h = 1
kin_visc_l = 1

mu_h = rho_h * kin_visc_h
mu_l = rho_l * kin_visc_l
mu_ratio = mu_l / mu_h

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

reference_lattice_size = 32
gauge = 2
lattice_size = int(gauge * reference_lattice_size)

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'IBB')

# ibb_Poiseuille_config_P00_00000000

case_folder = os.path.join(main_folder, f"Poiseuille_2_phase_v100")
oldest = find_oldest_iteration(case_folder)

# filename_vtk = f'{CollisionType}_ux_{ux}_k_{k}_sigma_{Sigma02}_size_{lattice_size}lu_VTK_P00_{oldest}.vti'
# filename_vtk = f'ibb_Poiseuille_VTK_P00_{oldest}.vti'
# filename_vtk = f'cm_Poisseulle_benchmark_VTK_P00_{oldest}.vti'

filename_vtk = get_vti_from_iteration(case_folder, oldest)

filepath_vtk = os.path.join(case_folder, filename_vtk)
vti_reader = VTIFile(filepath_vtk)
# T_num = vti_reader.get("T")
# T_num_slice = T_num[:, :, 1]

(U_num_x, _, _) = vti_reader.get("U", is_vector=True)
# U_num_x_slice = U_num_x[:, int(reference_lattice_size/2), 1]
U_num_x_slice = U_num_x[:, int(reference_lattice_size/2)]  # 2d?

# ySIZE = lattice_size
# xSIZE = lattice_size

ySIZE = U_num_x_slice.shape[0]


y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5

# -------- anal solution ---------------

h = 49  # distance from the center to the channel walls
y_anal = np.linspace(-h, h, 2*h, endpoint=False) + 0.5
uc = 0.001

# gx = 4.04e-06
gx = calc_gx(uc, mu_l, mu_h, rho_l, rho_h, h)
poiseuilleAnal = TwoPhasePoiseuilleAnal(gx=gx, mu_l=mu_l, mu_h=mu_h, rho_h=rho_h, rho_l=rho_l, h=h)
u_anal = np.array([poiseuilleAnal.get_u_profile(y_anal[i]) for i in range(len(y_anal))])


y_fd = np.linspace(-h, h, 10000, endpoint=False)
poiseuille_fd = TwoPhasePoiseuilleFD(gx=gx, mu_l=mu_l, mu_h=mu_h, rho_h=rho_h, rho_l=rho_l, h=h)
u_fd = poiseuille_fd.get_u_profile(y_fd, W=5)
###################################################################################################################

fig_name = f'Poiseuille_anal_vs_fd_rho{rho_h/rho_l}_v{kin_visc_h/kin_visc_l}.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

# plt.plot(u_anal, y_anal+len(y_anal)/2+1,
plt.plot(u_anal, y_anal + len(y_grid) / 2,
         color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
         label=r'$analytical \, solution$')

plt.plot(u_fd, y_fd + len(y_grid) / 2,
         color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
         label=r'$FD \, solution$')

plt.plot(U_num_x_slice, y_grid,
# plt.plot(U_num_x_slice, y_grid,
# plt.plot(U_num_x_slice, y_grid-len(y_grid)/2,
         color="black", marker="x", markevery=1, markersize=5, linestyle="", linewidth=2,
         label=r'$LBM \, solution$')

# ------ format y axis ------ #
yll = y_grid.min()
yhl = y_grid.max()
# axes.set_ylim([yll, yhl])
# axes.set_yticks(np.linspace(yll, yhl, 8))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.set_yticks([0.5, 1.5, 2.5, 31.5, 32, 32.5, 61.5, 62.5, 63.5])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
# ------ format x axis ------ #
# plt.xlim(0, int(xSIZE / 2))
# plt.xlim(int(xSIZE / 2), xSIZE)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'two phase Poiseuille flow')

plt.xlabel(r'$u_x$')
plt.ylabel(r'$y$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure
