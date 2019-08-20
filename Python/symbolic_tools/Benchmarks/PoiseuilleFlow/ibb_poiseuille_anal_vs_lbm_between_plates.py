import matplotlib.pyplot as plt
import numpy as np
import os
from Benchmarks.PoiseuilleFlow.PoiseuilleAnal import calc_gx_between_plates, OnePhasePoiseuilleAnalBetweenPlates
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker
import pwd
from DataIO.VTIFile import VTIFile
import matplotlib.pylab as pylab

rho = 1
kin_visc = 0.1
mu = rho * kin_visc

diam = 15
q = 0.75
expected_wall_location = 1.5 - q
effdiam = diam - 2*expected_wall_location
uc = 0.01

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_(I)BB_Poiseuille')


def read_ux(folder):
    folder = strip_folder_name(folder)
    case_folder = os.path.join(main_folder, folder)
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest)
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk)

    (u_num_x, _, _) = vti_reader.get("U", is_vector=True)
    ux_slice = u_num_x[:, 1, 1]
    return ux_slice


U_ibb_num_x_slice = read_ux(f"ibb_Poiseuille_nu_{kin_visc}_effdiam_{effdiam}")
U_bb_num_x_slice = read_ux(f"bb_Poiseuille_nu_{kin_visc}_effdiam_{effdiam}")
# U_bb_num_x_slice = read_ux(f"bb_Poiseuille_nu_{kin_visc}_effdiam_{diam - 2}")

ySIZE = U_ibb_num_x_slice.shape[0]
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5


# -------- anal solution ---------------
y_anal = np.arange(q, effdiam, 1)
y_anal = np.concatenate(([0], y_anal, [effdiam]))
gx = calc_gx_between_plates(uc, mu, mu, rho, rho, effdiam / 2)
poiseuilleAnal = OnePhasePoiseuilleAnalBetweenPlates(gx=gx, nu=kin_visc, H=effdiam)
u_anal = np.array([poiseuilleAnal.get_u_profile(y_anal[i]) for i in range(len(y_anal))])


###################################################################################################################
if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/Poiseuille_anal_vs_lbm_v{eat_dots_for_texmaker(kin_visc)}_q{eat_dots_for_texmaker(q)}_effdiam_{eat_dots_for_texmaker(effdiam)}.pdf'

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (14, 8),
          'axes.labelsize': 'xx-large',
          'axes.titlesize': 'xx-large',
          'xtick.labelsize': 'xx-large',
          'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

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

plt.plot(U_ibb_num_x_slice, y_grid,
         color="black", marker="x", markevery=1, markersize=7, linestyle="", linewidth=2,
         label=r'$LBM \, IBB$')

plt.plot(U_bb_num_x_slice, y_grid,
         color="black", marker="v", markevery=1, markersize=6, linestyle="", linewidth=2,
         label=r'$LBM \, BB$')

# ------ format y axis ------ #
yll = y_grid.min()
yhl = y_grid.max()
axes.set_ylim([0, 15])
# axes.set_yticks(np.linspace(yll, yhl, 8))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.set_yticks([0.5, 1.5, 2.5, 31.5, 32, 32.5, 61.5, 62.5, 63.5])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
# ------ format x axis ------ #
plt.xlim(0, 0.011)
# plt.xlim(int(xSIZE / 2), xSIZE)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

title = f'IBB Poiseuille flow:\n' + r'$ \nu = $' + f'{kin_visc:.2e} ' + r'$D_{eff}=$' + f'{effdiam:.2e}, q={q:.2e}'
title = ''  # skip title for .tex
plt.title(title)


plt.xlabel(r'$u_x$')
plt.ylabel(r'$y$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
# plt.show()

plt.close(fig)  # close the figure
