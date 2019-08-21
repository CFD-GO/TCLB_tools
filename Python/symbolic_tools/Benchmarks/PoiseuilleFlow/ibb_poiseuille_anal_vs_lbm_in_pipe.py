import matplotlib.pyplot as plt
import numpy as np
import os
from Benchmarks.PoiseuilleFlow.PoiseuilleAnal import calc_gx_in_pipe, OnePhasePoiseuilleAnalInPipe
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, get_r_from_xy
from DataIO.helpers import calc_L2, calc_mse
import pwd
from DataIO.VTIFile import VTIFile
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.VTIFile import VTIFile


rho = 1.
kin_visc = 0.1
mu = rho * kin_visc

effdiam = 28
uc = 0.01

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IBB_PoiseuillePipe')

case_folder = f'ibb_PoiseuillePipe_CM_HIGHER_nu_{kin_visc}_effdiam_{effdiam}'


def read_data_from_LBM(case_folder):
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    T_num = vti_reader.get("T")
    [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    ny, nx, nz = T_num.shape

    uz_num_slice = uz_num[:, :, 1]
    # T_num_slice = T_num[:, :, 1]
    #
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return uz_num_slice


uz_num_slice = read_data_from_LBM(os.path.join(main_folder, case_folder))
ny, nx = uz_num_slice.shape

# -------- anal solution ---------------
# y_anal = np.arange(q, effdiam, 1)
# y_anal = np.concatenate(([0], y_anal, [effdiam]))
# r_anal = np.arange(0, effdiam/2, 1)
gx = calc_gx_in_pipe(uc, rho, kin_visc, D=effdiam)
poiseuilleAnal = OnePhasePoiseuilleAnalInPipe(gx=gx, nu=kin_visc, D=effdiam)
# u_anal = np.array([poiseuilleAnal.get_u_profile(r_anal[i]) for i in range(len(r_anal))])


x0 = 63.5
y0 = 63.5

x_grid = np.linspace(0, nx, ny, endpoint=False) + 0.5
y_grid = np.linspace(0, nx, ny, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)

u_anal = np.zeros((ny, nx))
r_anal = np.zeros((ny, nx))
cuttoff = effdiam / 2. - 1
for i in range(ny):
    for j in range(nx):
        r = get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
        r_anal[i, j] = r
        u_anal[i, j] = poiseuilleAnal.get_u_profile(r)
        if r > cuttoff:
            u_anal[i, j] = np.nan


not_nan_mask = ~np.isnan(u_anal)
u_anal_masked = u_anal[not_nan_mask]
uz_num_slice_masked = uz_num_slice[not_nan_mask]
r_anal_masked = r_anal[not_nan_mask]

uz_mse = calc_mse(u_anal_masked, uz_num_slice_masked)
uz_L2 = calc_L2(u_anal_masked, uz_num_slice_masked)

# print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
# print(f"T_L2={T_L2[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")


###################################################################################################################
def cntr_plot():
    if not os.path.exists('plots'):
        os.makedirs('plots')

    fig_name = f'plots/Poiseuille_pipe_anal_vs_lbm_v{eat_dots_for_texmaker(kin_visc)}_effdiam_{eat_dots_for_texmaker(effdiam)}.png'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)


    fig = plt.figure(figsize=(12, 8))
    # ax = fig.gca(projection='3d')
    ax = fig.gca()
    # alpha=1, rstride=1, cstride=1)
    # ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5, label='T_err_field')
    # ax.plot_surface(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, label='T_num')
    # ax.plot_surface(xx, yy, T_anal,  cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.1, label='T_anal')
    # ax.contourf(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, label='T_num')

    cntr = ax.pcolormesh(xx, yy, uz_num_slice - u_anal, cmap='coolwarm', label='U_num')  # this one has smooth colors
    # cntr = ax.contourf(xx, yy, uz_num_slice, cmap='coolwarm', antialiased=True)  # this one is has step colors

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    ax.set_aspect('equal')

    # Customize the z axis.
    # ax.set_zlim(-.1, 1.05)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter(FormatStrFormatter('%.1e'))


    # Add a color bar which maps values to colors.

    fig.colorbar(cntr, shrink=0.5, aspect=5)

    title = f'IBB Poiseuille flow:\n' + r'$ \nu = $' + f'{kin_visc:.2e} ' + r'$D_{eff}=$' + f'{effdiam:.2e}'
    title = ''  # skip title for .tex
    plt.title(title)

    # Major ticks every 20, minor ticks every 5
    major_ticks = np.arange(0, nx, 10)
    minor_ticks = np.arange(0, nx, 1)

    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)

    # And a corresponding grid
    ax.grid(which='both')

    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    # plt.grid(True)

    # fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
    # ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)

    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/Plot_from_3D_data.png'

    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()
    # plt.close(fig)  # close the figure


def slice_plot():
    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/Poiseuille_pipe_anal_vs_lbm_v{eat_dots_for_texmaker(kin_visc)}_effdiam_{eat_dots_for_texmaker(effdiam)}.png'

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

    plt.plot(u_anal, r_anal,
             color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
             label=r'$analytical \, solution$')

    # plt.plot(u_anal, y_anal + len(y_grid) / 2,
    #          color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
    #          label=r'$analytical \, solution$')

    # plt.plot(u_fd, y_fd + len(y_grid) / 2,
    #          color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
    #          label=r'$FD \, solution$')

    plt.plot(uz_num_slice, r_anal,
             color="black", marker="x", markevery=1, markersize=7, linestyle="", linewidth=2,
             label=r'$LBM \, IBB$')
    #
    # plt.plot(U_bb_num_x_slice, y_grid,
    #          color="black", marker="v", markevery=1, markersize=6, linestyle="", linewidth=2,
    #          label=r'$LBM \, BB$')
    #
    # # ------ format y axis ------ #
    # yll = y_grid.min()
    # yhl = y_grid.max()
    axes.set_ylim([0, 2. + effdiam/2.])
    # axes.set_yticks(np.linspace(yll, yhl, 8))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
    # axes.set_yticks([0.5, 1.5, 2.5, 31.5, 32, 32.5, 61.5, 62.5, 63.5])
    # axes.yaxis.set_major_formatter(xfmt)

    # plt.yscale('log')
    # ------ format x axis ------ #
    # plt.xlim(0, 0.011)
    # plt.xlim(int(xSIZE / 2), xSIZE)

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    title = f'IBB Poiseuille flow:\n' + r'$ \nu = $' + f'{kin_visc:.2e} ' + r'$D_{eff}=$' + f'{effdiam:.2e}'
    title = ''  # skip title for .tex
    plt.title(title)


    plt.xlabel(r'$u_x$')
    plt.ylabel(r'$r$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    # plt.show()

    # plt.close(fig)  # close the figure

cntr_plot()
# # 2D clip
# r_anal = r_anal[:, 63]
# u_anal = u_anal[:, 63]
# uz_num_slice = uz_num_slice[:, int(nx / 2)]
# slice_plot()

print('bye')
