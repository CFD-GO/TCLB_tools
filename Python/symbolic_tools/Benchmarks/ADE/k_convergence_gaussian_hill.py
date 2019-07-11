

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

C0 = 1
X0 = Matrix([lattice_size/2, lattice_size/2])
U = Matrix([0, 0])
conductivities = ['0.01', '0.001', '1e-04', '1e-05']
Sigmas = [50, 75, 100]


for Sigma02 in Sigmas:
    def get_t_err(main_folder, collision_type):
        n = len(conductivities)
        T_mse = np.zeros(n)
        T_L2 = np.zeros(n)

        time_spot = None

        for g in range(n):
            folder = os.path.join(main_folder,
                                       f"{collision_type}_ux_{U[0]}_k_{conductivities[g]}_sigma_{Sigma02}_size_{lattice_size}lu")

            oldest = find_oldest_iteration(folder)
            time_spot = int(oldest)
            filename_vtk = f'{collision_type}_ux_{U[0]}_k_{conductivities[g]}_sigma_{Sigma02}_size_{lattice_size}lu_VTK_P00_{oldest}.vti'

            filepath_vtk = os.path.join(folder, filename_vtk)
            vti_reader = VTIFile(filepath_vtk)
            T_num = vti_reader.get("T")
            T_num_slice = T_num[:, :, 1]

            ySIZE, xSIZE = T_num_slice.shape
            assert ySIZE == xSIZE == lattice_size

            def calc_anal_solution():
                ySIZE = lattice_size
                xSIZE = lattice_size

                x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
                y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
                xx, yy = np.meshgrid(x_grid, y_grid)

                total_time = 1

                T_anal = np.zeros((ySIZE, xSIZE, total_time))

                gha = GaussianHillAnal2D(C0, X0, U, Sigma02, float(conductivities[g]))
                for i in range(ySIZE):
                    print(f"running i/ySIZE = {i}/{ySIZE}...")
                    for j in range(xSIZE):
                        T_anal[i][j][0] = gha.get_concentration(Matrix([xx[i][j], yy[i][j]]), time_spot)  # lets cheat

                T_anal = T_anal[:, :, 0]  # take time slice
                return xx, yy, T_anal


            shall_recalculate_results = False
            dump_fname = os.path.join(main_folder, f'dumps',
                                      f'ux_{U[0]}_k_{conductivities[g]}_sigma_{Sigma02}_size_{lattice_size}_time_spot_{oldest}.npy')

            if os.path.isfile(dump_fname) and not shall_recalculate_results:
                print(f'{dump_fname} found, loading results from disc')
                (xx, yy, T_anal) = np.load(dump_fname)
            else:
                print(f'recalculating results for: {dump_fname}')
                xx, yy, T_anal = calc_anal_solution()
                dump_folder = os.path.dirname(dump_fname)
                if not os.path.exists(dump_folder):
                    os.makedirs(dump_folder)
                np.save(dump_fname, (xx, yy, T_anal))

            T_err_field = T_anal - T_num_slice
            T_L2[g] = calc_L2(T_anal, T_num_slice)
            T_mse[g] = calc_mse(T_anal, T_num_slice)

            # print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
            print(f"T_L2={T_L2[g]:.2e} for k = {conductivities[g]}")

        return T_L2, time_spot
        # return T_mse


    home = pwd.getpwuid(os.getuid()).pw_dir
    # T_err_CM = get_t_err(home, 'CM')

    T_err_CM_HIGHER, time_spot = get_t_err(os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_advection'), 'CM_HIGHER')
    T_err_CM, _ = get_t_err(os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianHill_advection'), 'CM')

    fconductivities = [float(c) for c in conductivities]
    print("------------------------------------ PLOT ------------------------------------")
    fig_name = f'ADE_GaussianHill_lattice={lattice_size}[lu]_sig={Sigma02}_Ux={U[0]}_time={time_spot}.png'
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(14, 8))

    axes = plt.gca()

    plt.plot(fconductivities, T_err_CM_HIGHER,
             color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
             label=f'CM HIGHER')

    plt.plot(fconductivities, T_err_CM,
             color="black", marker="x", markevery=1, markersize=7, linestyle="", linewidth=2,
             label=f'CM')


    # ------ format y axis ------ #
    # yll = T_err_CM_HIGHER.min()
    # yhl = T_err_CM_HIGHER.max()
    # axes.set_ylim([0, yhl])

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

    plt.xlabel(r'$conductivity$', fontsize=18)
    plt.ylabel(r'$T: \; L_2 \, error \, norm $', fontsize=18)
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure
