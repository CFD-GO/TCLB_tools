import matplotlib.pylab as pylab
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.pyplot as plt
import os
import numpy as np


def cntr_plot(anal_field, num_field, xx, yy, conductivity, eff_pipe_diam):
    print("---------- PLOTTING -------------")
    fig = plt.figure(figsize=(12, 8))
    # ax = fig.gca(projection='3d')
    ax = fig.gca()
    # alpha=1, rstride=1, cstride=1)
    # ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5, label='T_err_field')
    # ax.plot_surface(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, label='T_num')
    # ax.plot_surface(xx, yy, T_anal,  cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.1, label='T_anal')
    # ax.contourf(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, label='T_num')

    cntr = ax.pcolormesh(xx, yy, num_field - anal_field, cmap='coolwarm', label='T_num')  # this one has smooth colors
    # cntr = ax.contourf(xx, yy, T_num_slice, cmap='coolwarm', antialiased=True)  # this one is has step colors

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

    title = f'iabb_ruraWrurze\n' \
            f'_conducivity={conductivity} eff_pipe_diam={eff_pipe_diam}'
    plt.title(title)

    # # Major ticks every 20, minor ticks every 5
    # major_ticks = np.arange(0, nx, nx/5)
    # minor_ticks = np.arange(0, nx, nx/10)
    #
    # ax.set_xticks(major_ticks)
    # ax.set_xticks(minor_ticks, minor=True)
    # # ax.set_yticks(major_ticks)
    # # ax.set_yticks(minor_ticks, minor=True)
    #
    # # And a corresponding grid
    # ax.grid(which='both')
    #
    # # Or if you want different settings for the grids:
    # ax.grid(which='minor', alpha=0.2)
    # ax.grid(which='major', alpha=0.5)
    plt.grid(True)  # or use default grid

    # fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
    # ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)

    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/{title}.png'

    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()
    plt.close(fig)  # close the figure
    # print('bye')