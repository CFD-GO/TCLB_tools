from Benchmarks.ADE.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeNeumann
from DataIO.helpers import find_oldest_iteration, calc_mse, calc_L2
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pwd
from DataIO.VTIFile import VTIFile
from Benchmarks.ADE.GaussianHillAnal2D import GaussianHillAnal2D

from sympy.matrices import Matrix, diag
# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

reference_lattice_size = 32
gauge = 8
lattice_size = int(gauge * reference_lattice_size)


##### prepare anal solultion ####

C0 = 1
X0 = Matrix([lattice_size/2, lattice_size/2])
U = Matrix([0, 0])
Sigma02 = 100
k = '0.001'
# conductivities = ['0.01', '0.001', '1e-04', '1e-05']

gha = GaussianHillAnal2D(C0, X0, U, Sigma02, float(k))

home = pwd.getpwuid(os.getuid()).pw_dir

CollisionType ="CM_HIGHER"


main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_advection',
                           f"{CollisionType}_ux_{U[0]}_k_{k}_sigma_{Sigma02}_size_{lattice_size}lu")

oldest = find_oldest_iteration(main_folder)
# oldest = "00001000"
filename_vtk = f'{CollisionType}_ux_{U[0]}_k_{k}_sigma_{Sigma02}_size_{lattice_size}lu_VTK_P00_{oldest}.vti'

filepath_vtk = os.path.join(main_folder, filename_vtk)
vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")
T_num_slice = T_num[:, :, 1]


ySIZE = lattice_size
xSIZE = lattice_size

x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)

total_time = 1
time_spot = int(oldest)
T_anal = np.zeros((ySIZE, xSIZE, total_time))


for i in range(ySIZE):
    print(f"running i/ySIZE = {i}/{ySIZE}...")
    for j in range(xSIZE):

        T_anal[i][j][0] = gha.get_concentration(Matrix([xx[i][j], yy[i][j]]), time_spot)  # lets cheat
        # for t in range(total_time):
            # T_anal[i][j][t] = gha.get_concentration(Matrix([xx[i][j], yy[i][j]]), t)


T_anal = T_anal[:, :, 0]  # take time slice

T_err_field = T_anal - T_num_slice
T_L2 = calc_L2(T_anal, T_num_slice)
T_mse = calc_mse(T_anal, T_num_slice)

###################################################################################################################
shall_clip_to_2D = True

if shall_clip_to_2D:
    T_anal = T_anal[:, int(ySIZE / 2)]  # half X slice
    T_num_slice = T_num_slice[:, int(ySIZE / 2)]  # half X slice
    T_err_field = T_err_field[:, int(ySIZE / 2)]
    # -------------------- make dummy plot --------------------
    fig_name = f'ADE_GaussianHill_lattice={lattice_size}[lu]_sig={Sigma02}_Ux={U[0]}_k={k}_time={time_spot}_2Dslice.png'
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(14, 8))

    axes = plt.gca()

    plt.plot(x_grid, T_anal,
             color="black", marker="", markevery=5, markersize=5, linestyle="-", linewidth=2,
             label=r'$analytical \, solution$')

    plt.plot(x_grid, T_num_slice,
             color="black", marker="", markevery=5, markersize=5, linestyle="-.", linewidth=2,
             label=r'$numerical \, solution$')

    plt.plot(x_grid, T_err_field,
             color="black", marker="", markevery=5, markersize=5, linestyle="-.", linewidth=2,
             label=r'$error \, field$')
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
    # ax.plot_surface(xx, yy, T_anal, cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.01, label='T_anal')
    # ax.plot_surface(xx, yy, T_num_slice, cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, alpha=1, label='T_num')
    ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.25, alpha=1, label='T_err_field')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('T')

    # Customize the z axis.
    # ax.set_zlim(-.1, .1)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.title(f'Advection - Diffusion of Gaussian Hill at t={time_spot}'
              f'\n' r'$T_{L2}$=' + f'{T_L2:.2e}')

    plt.grid(True)

    fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
    ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
    plt.show()

    print("bye")

