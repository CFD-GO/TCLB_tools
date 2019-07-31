import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import os
from Benchmarks.TwoPhasePoiseuilleFlow.TwoPhasePoiseuilleAnal import TwoPhasePoiseuilleAnal, calc_gx, OnePhasePoiseuilleAnal
from Benchmarks.TwoPhasePoiseuilleFlow.TwoPhasePoiseuilleFD import TwoPhasePoiseuilleFD
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, calc_mse, calc_L2
import pwd
from DataIO.VTIFile import VTIFile


rho_h = 1
rho_l = 1

kin_visc_h = 0.1
kin_visc_l = kin_visc_h

mu_h = rho_h * kin_visc_h
mu_l = rho_l * kin_visc_l
mu_ratio = mu_l / mu_h

diam = 63
q = 0.75
expected_wall_location = 1.5 - q
effdiam = diam - 2*expected_wall_location
uc = 0.01

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IBB_Poiseuille')

def read_ux(folder):
    case_folder = os.path.join(main_folder, folder)
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest)
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk)

    (u_num_x, _, _) = vti_reader.get("U", is_vector=True)
    ux_slice = u_num_x[:, 1, 1]
    return ux_slice


U_num_x_slice = read_ux(f"ibb_Poiseuille_nu_{kin_visc_h}_effdiam_{effdiam}")
ySIZE = U_num_x_slice.shape[0]
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5

U_bb_num_x_slice = read_ux(f"bb_Poiseuille_nu_{kin_visc_h}_effdiam_{diam-2.0}")

# -------- anal solution ---------------

# gx = calc_gx(uc, mu_l, mu_h, rho_l, rho_h, effdiam/2)
# y_anal = np.linspace(-31, 31, 62, endpoint=False)
# poiseuilleAnal = TwoPhasePoiseuilleAnal(gx=gx, mu_l=mu_l, mu_h=mu_h, rho_h=rho_h, rho_l=rho_l, r=effdiam / 2)
# u_anal = np.array([poiseuilleAnal.get_u_profile(y_anal[i]) for i in range(len(y_anal))])


y_anal = np.arange(q, effdiam, 1)
y_anal = np.concatenate(([0], y_anal, [effdiam]))
gx = calc_gx(uc, mu_l, mu_h, rho_l, rho_h, effdiam/2)
poiseuilleAnal = OnePhasePoiseuilleAnal(gx=gx, nu=kin_visc_h, D=effdiam)
u_anal = np.array([poiseuilleAnal.get_u_profile(y_anal[i]) for i in range(len(y_anal))])

# y_fd = np.linspace(-31, 31, 1000, endpoint=False)
# poiseuille_fd = TwoPhasePoiseuilleFD(gx=gx, mu_l=mu_l, mu_h=mu_h, rho_h=rho_h, rho_l=rho_l, r=effdiam / 2)
# u_fd = poiseuille_fd.get_u_profile(y_fd, W=5)
###################################################################################################################

fig_name = f'Poiseuille_anal_vs_lbm_v{kin_visc_h/kin_visc_l}_q{q}_effdiam_{effdiam}.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

plt.plot(u_anal, y_anal + expected_wall_location,
         color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
         label=r'$analytical \, solution$')

# plt.plot(u_anal, y_anal + len(y_grid) / 2,
#          color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
#          label=r'$analytical \, solution$')

# plt.plot(u_fd, y_fd + len(y_grid) / 2,
#          color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
#          label=r'$FD \, solution$')

plt.plot(U_num_x_slice, y_grid,
         color="black", marker="x", markevery=1, markersize=7, linestyle="", linewidth=2,
         label=r'$LBM \, IBB$')

plt.plot(U_bb_num_x_slice, y_grid,
         color="black", marker="v", markevery=1, markersize=6, linestyle="", linewidth=2,
         label=r'$LBM \, BB$')

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


plt.title(f'IBB Poiseuille flow:\n' +
          r'$ \nu = $' + f'{kin_visc_h/kin_visc_l} ' + r'$D_{eff}=$' + f'{effdiam}, q={q}')

plt.xlabel(r'$u_x$')
plt.ylabel(r'$y$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure
