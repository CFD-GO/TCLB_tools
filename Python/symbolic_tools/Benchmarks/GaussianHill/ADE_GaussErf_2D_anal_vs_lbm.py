import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import os
from Benchmarks.GaussianHill.GaussianHillAnal2D import GaussianHillAnal2D

from sympy.matrices import Matrix

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

reference_lattice_size = 32
gauge = 2
lattice_size = int(gauge * reference_lattice_size)


##### prepare anal solultion ####

C0 = 1
X0 = Matrix([lattice_size/2, lattice_size/2])
U = Matrix([0.1, 0])
Sigma02 = 4
k = 0.166666
gha = GaussianHillAnal2D(C0, X0, U, Sigma02, k)



ySIZE = lattice_size
xSIZE = lattice_size

x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)

total_time = 1
time_spot = 0
T_anal = np.zeros((ySIZE, xSIZE, total_time))


for i in range(ySIZE):
    for j in range(xSIZE):
        T_anal[i][j][0] = gha.get_concentration(Matrix([xx[i][j], yy[i][j]]), time_spot)  # lets cheat
        # for t in range(total_time):
            # T_anal[i][j][t] = gha.get_concentration(Matrix([xx[i][j], yy[i][j]]), t)


T_anal = T_anal[:, :, time_spot]  # take time slice
###################################################################################################################
shall_clip_to_2D = True

if shall_clip_to_2D:

    T_anal = T_anal[:, int(ySIZE / 2)]  # half X slice



    fig_name = f'ADE_GaussianHill_lattice={lattice_size}[lu]_sig={Sigma02}_Ux={U[0]}_k={k}_time={time_spot}_2D.png'

    # -------------------- make dummy plot --------------------
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(14, 8))

    axes = plt.gca()

    plt.plot(x_grid, T_anal,
             color="black", marker="", markevery=5, markersize=5, linestyle="-", linewidth=2,
             label=r'$analytical \, solution$')

    # plt.plot(u_fd, y_,
    #          color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
    #          label=r'$FD \, solution$')


    # ------ format y axis ------ #
    yll = T_anal.min()
    yhl = T_anal.max()
    axes.set_ylim([0, C0])
    # axes.set_yticks(np.linspace(yll, yhl, 5))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
    # axes.yaxis.set_major_formatter(xfmt)

    # plt.yscale('log')
    # ------ format x axis ------ #
    # plt.xlim(0, int(xSIZE / 2))
    # plt.xlim(int(xSIZE / 2), xSIZE)

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


    plt.title(f'Advection - Diffusion of Gaussian Hill at t={time_spot}')

    plt.xlabel(r'$x$')
    plt.ylabel(r'$Temperature$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    # plt.close(fig)  # close the figure

else:
    print("---------- PLOTTING -------------")

    fig_name = f'ADE_GaussianHill_lattice={lattice_size}[lu]_sig={Sigma02}_Ux={U[0]}_k={k}_time={time_spot}_3D.png'

    fig = plt.figure(figsize=(12, 8))
    ax = fig.gca(projection='3d')

    # alpha=1, rstride=1, cstride=1)
    # ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5,
    #                 label='T_err_field')  # coolwarm
    ax.plot_surface(xx, yy, T_anal, cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.01, label='T_anal')
    # ax.plot_surface(xx, yy, T_num, cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, alpha=1, label='T_num')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Temperature')

    # Customize the z axis.
    # ax.set_zlim(-.1, .1)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.title(f'Advection - Diffusion of Gaussian Hill at t={time_spot}')

    plt.grid(True)

    fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
    ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
    plt.show()

    print("bye")

