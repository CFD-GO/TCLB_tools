from DataIO.helpers import find_oldest_iteration, calc_mse, calc_L2
import matplotlib.pyplot as plt
import os
import pwd
from DataIO.VTIFile import VTIFile
from Benchmarks.GaussianHill.GaussianHillAnal2D import GaussianHillAnal2D, prepare_anal_data_ADE_Gaussion_Hill

from sympy.matrices import Matrix

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

reference_lattice_size = 32
gauge = 8
lattice_size = int(gauge * reference_lattice_size)

C0 = 1
X0 = Matrix([lattice_size/2, lattice_size/2])
ux = 0
# ux_str = str(Ux).rstrip('0')
Sigma02 = 100
conductivities = ['0.01', '0.001', '1e-04', '1e-05']

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_advection')

shall_clip_to_2D = False

for CollisionType in ["CM", "CM_HIGHER"]:
    for k in conductivities:
        case_folder = os.path.join(main_folder, f"{CollisionType}_ux_{ux}_k_{k}_sigma_{Sigma02}_size_{lattice_size}lu")

        oldest = find_oldest_iteration(case_folder)
        filename_vtk = f'{CollisionType}_ux_{ux}_k_{k}_sigma_{Sigma02}_size_{lattice_size}lu_VTK_P00_{oldest}.vti'
        filepath_vtk = os.path.join(case_folder, filename_vtk)
        vti_reader = VTIFile(filepath_vtk)
        T_num = vti_reader.get("T")
        T_num_slice = T_num[:, :, 1]

        ySIZE = lattice_size
        xSIZE = lattice_size

        dump_file_path = os.path.join(main_folder, f'dumps',
                                      f'ux_{ux}_k_{k}_sigma_{Sigma02}_size_{lattice_size}_time_spot_{oldest}.npy')

        gha = GaussianHillAnal2D(C0, X0, Sigma02, float(k))
        xx, yy, T_anal = prepare_anal_data_ADE_Gaussion_Hill(gha, ux, oldest, ySIZE, xSIZE, dump_file_path,
                                                             shall_recalculate_results=False)

        T_err_field = T_anal - T_num_slice
        T_L2 = calc_L2(T_anal, T_num_slice)
        T_mse = calc_mse(T_anal, T_num_slice)

        ###################################################################################################################

        if shall_clip_to_2D:
            T_anal = T_anal[int(ySIZE / 2), :]  # half X slice
            T_num_slice = T_num_slice[int(ySIZE / 2), :]  # half X slice
            T_err_field = T_err_field[int(ySIZE / 2), :]

            x_grid = xx[0, :]
            # -------------------- make dummy plot --------------------
            fig_name = f'ADE_GaussianHill_lattice={lattice_size}[lu]_sig={Sigma02}_ux={ux}_k={k}_time={int(oldest)}_2Dslice.png'
            plt.rcParams.update({'font.size': 14})
            plt.figure(figsize=(14, 8))

            axes = plt.gca()

            plt.plot(x_grid, T_anal,
                     color="black", marker="", markevery=5, markersize=5, linestyle="-", linewidth=2,
                     label=r'$analytical \, solution$')

            plt.plot(x_grid, T_num_slice,
                     color="black", marker="x", markevery=5, markersize=5, linestyle="-.", linewidth=2,
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

            plt.title(f'Advection - Diffusion of Gaussian Hill at t={int(oldest)}'
                      f'\n' r'$T_{L2}$=' + f'{T_L2:.2e}')

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

            fig_name = f'ADE_GaussianHill_{CollisionType}_ux={ux}_sig={Sigma02}_k={float(k):.0e}_time={int(oldest)}_lattice={lattice_size}[lu]_contour.png'
            plt.rcParams.update({'font.size': 24})
            fig = plt.figure(figsize=(12, 8))


            # ax = fig.gca(projection='3d')
            # alpha=1, rstride=1, cstride=1)
            # ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5,
            #                 label='T_err_field')  # coolwarm

            # ax.plot_surface(xx, yy, 0.1+T_anal, cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.01, label='T_anal')
            # ax.plot_surface(xx, yy, T_num_slice, cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, alpha=1, label='T_num')
            # ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.25, alpha=1, label='T_err_field')
            # Add a color bar which maps values to colors.
            # fig.colorbar(surf, shrink=0.5, aspect=5)

            ax = fig.gca()
            # cntr = ax.pcolormesh(xx, yy, T_num_slice, cmap='coolwarm', label='T_num')  # this one has smooth colors
            cntr = ax.pcolormesh(xx, yy, T_err_field, cmap='coolwarm', label='T_num')  # this one has smooth colors
            # cntr = ax.contourf(xx, yy, T_num_slice, cmap='coolwarm', antialiased=True)  # this one is has step colors
            # fig.colorbar(cntr, shrink=0.5, aspect=5)
            fig.colorbar(cntr, format='%.0e').set_label(label='', size=20)

            ax.set_xlabel('X', fontsize=24)
            ax.set_ylabel('Y', fontsize=24)
            # ax.set_zlabel('T')

            # Customize the z axis.
            # ax.set_zlim(-.1, .1)
            # ax.zaxis.set_major_locator(LinearLocator(10))

            # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))



            # plt.title(f'Advection - Diffusion of Gaussian Hill at t={int(oldest)}'
            #           f'\n' r'$T_{L2}$=' + f'{T_L2:.2e}')

            # plt.grid(True)

            # fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
            # ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
            fig = plt.gcf()  # get current figure
            fig.savefig(fig_name, bbox_inches='tight')
            plt.show()

            plt.close(fig)  # close the figure

            print("bye")

