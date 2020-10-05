from DataIO.helpers import find_oldest_iteration, calc_mse, calc_L2
import matplotlib.pyplot as plt
import numpy as np
import os
import pwd
import re
from fractions import Fraction
from DataIO.VTIFile import VTIFile
from DataIO.helpers import peel_the_skin_v2
from Benchmarks.GaussianHill.GaussianHillAnal2D import GaussianHillAnal2D, prepare_anal_data_ADE_Gaussion_Hill

from sympy.matrices import Matrix

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir

# time_SI = 10
# conductivity_SI = 2.0
# iterations = 2*10*np.array([6, 10, 100, 1000, 10000])
# main_dir = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_acoustic_scalling', '120-200000_iterations')

time_SI = 50
conductivity_SI = 7.0
iterations = 7*50*np.array([6, 10, 100, 1000, 10000])
main_dir = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_acoustic_scalling', '2100-3500000_iterations')

# time_SI = 100
# conductivity_SI = 4.0
# iterations = 4*100*np.array([6, 10, 100, 1000, 10000])
# main_dir = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_acoustic_scalling', '2400-4000000_iterations')

domain_size_SI = 256.0
lattice_size = 256

dx = domain_size_SI/lattice_size
conductivities = np.array((time_SI * conductivity_SI) / (iterations * dx * dx))
str_conductivities = [re.sub(r"/", 'o', str(Fraction(conductivity).limit_denominator(max_denominator=100000))) for conductivity in conductivities]

# Gaussian Hill coordinates
C0 = 1.
X0 = Matrix([lattice_size/2., lattice_size/2.])

# for ux in [0, 0.1]:
#     for Sigma02 in [50, 75, 100]:

plot_dir = 't_convergence_plots'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)


def plot_err_field(T_err_field, xx, yy, fig_name):
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
    cntr = ax.pcolormesh(xx, yy, T_err_field, cmap='coolwarm', label='T_num', shading='auto')  # this one has smooth colors
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


def plot_t_convergence_acoustic_scalling(conductivities, T_err_L2_BGK, T_err_L2_CM, T_err_L2_CM_HIGHER, T_err_L2_Cumulants, fig_name):
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(14, 8))
    # '{:.0e}'.format(x)
    axes = plt.gca()

    plt.plot(conductivities, T_err_L2_BGK,
             color="black", marker="^", markevery=1, markersize=7, linestyle="", linewidth=2,
             label=f'BGK')

    plt.plot(conductivities, T_err_L2_CM,
             color="black", marker="v", markevery=1, markersize=7, linestyle="", linewidth=2,
             label=f'CM')

    plt.plot(conductivities, T_err_L2_CM_HIGHER,
             color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
             label=f'CM TRT')

    plt.plot(conductivities, T_err_L2_Cumulants,
             color="black", marker="x", markevery=1, markersize=7, linestyle="", linewidth=2,
             label=f'Cumulants')

    # ------ format y axis ------ #
    # yll = T_err_CM_HIGHER.min()
    # yhl = T_err_CM_HIGHER.max()

    # axes.set_ylim([0.0004, 0.1])

    # axes.set_yticks(np.linspace(yll, yhl, 5))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
    # axes.yaxis.set_major_formatter(xfmt)

    plt.yscale('log')
    # ------ format x axis ------ #
    plt.xscale('log')
    # plt.xlim(0, int(xSIZE / 2))
    # plt.xlim(int(xSIZE / 2), xSIZE)

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    # plt.title(f'Advection - Diffusion of a Gaussian Hill at t={time_spot}')


    # normalnie
    # dk = np.logspace(np.log10(conductivities[0]), np.log10(conductivities[-1]), 100)
    # y = dk ** 1
    # y = y / y[0] * T_err_L2_BGK[0]
    # plt.loglog(dk, y, label=r'${x}$')
    #
    # y2 = dk ** 2
    # y2 = y2 / y2[0] * T_err_L2_BGK[0]
    # plt.loglog(dk, y2, label=r'${x^2}$')

    # odwrotnie
    # dk = np.logspace(np.log10(conductivities[-1]), np.log10(conductivities[0]), 100)
    # initial_error_1st = max(T_err_L2_CM)
    # y = min(dk) * initial_error_1st / dk
    # plt.loglog(dk, y, label=r'${-x}$')
    #
    # initial_error_2nd = max(T_err_L2_CM)
    # y2 = min(dk) * min(dk) * initial_error_2nd / (dk * dk)
    # plt.loglog(dk, y2, label=r'${-(x^2)}$')



    plt.xlabel(r'$k$', fontsize=22)
    plt.ylabel(r'$T: \; L_2 \, error \, norm $', fontsize=22)
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure


for ux in [0, 0.1]:
    for Sigma02 in [100]:
        def get_t_err(main_folder, collision_type):
            n = len(str_conductivities)
            T_mse = np.zeros(n)
            T_L2 = np.zeros(n)

            for g in range(n):
                folder = os.path.join(
                    main_folder,
                    f"{collision_type}_ux_{ux:.2e}_k_{str_conductivities[g]}_iterations_{iterations[g]}_sigma_{Sigma02}_size_{lattice_size}lu")

                oldest = find_oldest_iteration(folder)
                filename_vtk = f'{collision_type}_ux_{ux:.2e}_k_{str_conductivities[g]}_iterations_{iterations[g]}_sigma_{Sigma02}_size_{lattice_size}lu_VTK_P00_{oldest}.vti'

                filepath_vtk = os.path.join(folder, filename_vtk)
                vti_reader = VTIFile(filepath_vtk)
                T_num = vti_reader.get("T")
                T_num_slice = T_num[:, :, 1]

                ySIZE, xSIZE = T_num_slice.shape
                assert ySIZE == xSIZE == lattice_size

                dump_file_path = os.path.join(main_folder, f'dumps',
                                              f'ux_{ux}_k_{str_conductivities[g]}_iterations_{iterations[g]}_sigma_{Sigma02}_size_{lattice_size}_time_SI_{time_SI}.npy')


                assert lattice_size == domain_size_SI
                gha = GaussianHillAnal2D(C0, X0, Sigma02, conductivities[g])
                xx, yy, T_anal = prepare_anal_data_ADE_Gaussion_Hill(gha, ux, iterations[g], lattice_size, lattice_size, dump_file_path,
                                                                     shall_recalculate_results=False, reference_level=10.)
                T_err_field = T_anal - T_num_slice

                # alternatively
                # gha2 = GaussianHillAnal2D(C0, X0, Sigma02, conductivity_SI)
                # ux_SI = conductivity_SI/conductivities[g]*ux
                # xx, yy, T_anal2 = prepare_anal_data_ADE_Gaussion_Hill(gha2, ux_SI, time_SI, lattice_size, lattice_size,
                #                                                      dump_file_path,
                #                                                      shall_recalculate_results=True,
                #                                                      reference_level=10.)
                # T_err_field2 = T_anal2 - T_num_slice
                # xx = peel_the_skin_v2(xx, int(0.25 * lattice_size), int(0.75 * lattice_size))
                # yy = peel_the_skin_v2(yy, int(0.25 * lattice_size), int(0.75 * lattice_size))
                # T_anal = peel_the_skin_v2(T_anal, int(0.25*lattice_size), int(0.75*lattice_size))
                # T_num_slice = peel_the_skin_v2(T_num_slice, int(0.25 * lattice_size), int(0.75 * lattice_size))
                # T_err_field = T_anal - T_num_slice

                T_L2[g] = calc_L2(T_anal, T_num_slice)
                T_mse[g] = calc_mse(T_anal, T_num_slice)

                # print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
                print(f"{collision_type} T_L2={T_L2[g]:.5e} for k = {conductivities[g]}")

                print("------------------------------------ PLOT err field------------------------------------")
                fig_name = f'{plot_dir}/acoustic_scaling_GaussianHill_{collision_type}_ux={ux:.0e}_k_{str_conductivities[g]}_iterations_{iterations[g]}_sig={Sigma02}_time_SI={time_SI}_lattice={lattice_size}[lu]_err_field_contour.png'
                plot_err_field(T_err_field, xx, yy, fig_name)

            return T_L2
            # return T_mse



        T_err_L2_BGK = get_t_err(main_dir, 'BGK')
        T_err_L2_CM = get_t_err(main_dir, 'CM')
        T_err_L2_CM_HIGHER = get_t_err(main_dir, 'CM_HIGHER')
        T_err_L2_Cumulants = get_t_err(main_dir, 'Cumulants')
        print("------------------------------------ PLOT t convergence------------------------------------")
        fig_name = f'{plot_dir}/acoustic_scalling_GaussianHill_ux={ux:.0e}_sig={Sigma02}_time_SI={time_SI}_conductivity_SI_{conductivity_SI}_lattice={lattice_size}[lu]_t_convergence.png'
        plot_t_convergence_acoustic_scalling(conductivities, T_err_L2_BGK, T_err_L2_CM, T_err_L2_CM_HIGHER, T_err_L2_Cumulants, fig_name)
        # plot_t_convergence_acoustic_scalling(conductivities, T_err_L2_BGK, T_err_L2_BGK, T_err_L2_BGK, T_err_L2_BGK, fig_name)  # HACK:

print("Done.")
