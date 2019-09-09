import matplotlib.pyplot as plt
import numpy as np
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, get_r_from_xy
import matplotlib.pylab as pylab
import os
import matplotlib.ticker

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)

    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom

    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)


def cntr_plot(u_anal, uz_num_slice, kin_visc, effdiam, nx, xx, yy):
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

    cntr = ax.pcolormesh(xx, yy, u_anal, cmap='coolwarm', label='U_num')  # this one has smooth colors
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

    # fig.colorbar(cntr, shrink=0.5, aspect=5,  format=OOMFormatter(-6, mathText=False))

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


def slice_plot(u_anal, uz_num_slice, kin_visc, effdiam, r_anal):
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

