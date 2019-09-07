
import numpy as np
import os
from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import HeatConductionBetweenTwoPlates
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, calc_mse, calc_L2, strip_folder_name, eat_dots_for_texmaker
import pwd
from DataIO.VTIFile import VTIFile

from matplotlib.ticker import LinearLocator, FormatStrFormatter
from DataIO.VTIFile import VTIFile
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.pylab as pylab
rho = 1
conductivity = '0.01'
mu = rho * float(conductivity)

diameters = np.array([15, 31, 63])
qs = [0.25, 0.5, 0.75]
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
    T_num = vti_reader.get("T")
    T_num_slice = T_num[:, :, 1]
    # (u_num_x, _, _) = vti_reader.get("U", is_vector=True)
    # ux_slice = u_num_x[:, 1, 1]
    return T_num_slice


def delete_unphysical_data_from_wall_nodes(data):
    data = np.delete(data, 0, axis=0)
    data = np.delete(data, -1, axis=0)
    return data


for q in qs:
    n_diam = len(diameters)
    T_iabb_mse = np.zeros(n_diam)
    T_iabb_L2 = np.zeros(n_diam)
    T_abb_mse = np.zeros(n_diam)
    T_abb_L2 = np.zeros(n_diam)

    for d in range(n_diam):
        expected_wall_location = 1.5 - q
        effdiam = diameters[d] - 2 * expected_wall_location

        T_iabb_num_slice = read_ux(f"iabb_Poiseuille_k_{conductivity}_effdiam_{effdiam}")
        T_abb_num_slice = read_ux(f"abb_Poiseuille_k_{conductivity}_effdiam_{effdiam}")
        # ux_bb_num_slice = read_ux(f"bb_Poiseuille_nu_{kin_visc}_effdiam_{diameters[d] - 2}")

        ySIZE = T_iabb_num_slice.shape[0]
        y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5

        T_iabb_num_slice = delete_unphysical_data_from_wall_nodes(T_iabb_num_slice)
        T_abb_num_slice = delete_unphysical_data_from_wall_nodes(T_abb_num_slice)
        y_grid = delete_unphysical_data_from_wall_nodes(y_grid)

        # -------- anal solution ---------------
        y_anal = np.arange(q, effdiam, 1)
        # y_anal = np.concatenate(([0], y_anal, [effdiam])) # unphysical wall nodes
        hcbp = HeatConductionBetweenTwoPlates(T0=11, T2=10, Heff=effdiam)

        T_anal = np.array([hcbp.get_T_profile(y_anal[i]) for i in range(len(y_anal))])

        T_iabb_mse[d] = calc_mse(T_anal, T_iabb_num_slice)
        T_iabb_L2[d] = calc_L2(T_anal, T_iabb_num_slice)
        T_abb_mse[d] = calc_mse(T_anal, T_abb_num_slice)
        T_abb_L2[d] = calc_L2(T_anal, T_abb_num_slice)

        print(f"ux_mse={T_iabb_mse[d]:.2e} for k{conductivity}_q{q}_effdiam_{effdiam}")
        print(f"ux_L2={T_iabb_L2[d]:.2e} for k{conductivity}_q{q}_effdiam_{effdiam}")

    print("------------------------------------ Convergence  PLOT ------------------------------------")
    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/grid_convergence_conduction_between_plates_anal_vs_lbm_k{eat_dots_for_texmaker(conductivity)}_q{eat_dots_for_texmaker(q)}.png'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    initial_error_05st = 0.070
    y_05st = np.sqrt(diameters.min())*initial_error_05st/np.sqrt(diameters)

    # initial_error_1st = 0.18
    initial_error_1st = 1.15*max(np.concatenate((T_iabb_L2, T_abb_L2)))
    y_1st = diameters.min()*initial_error_1st/diameters
    # initial_error_2nd = 0.05
    initial_error_2nd = 0.85*min((max(T_iabb_L2), max(T_abb_L2)))
    y_2nd = diameters.min()*diameters.min()*initial_error_2nd/(diameters*diameters)

    fig1, ax1 = plt.subplots(figsize=(14, 8))
    plt.rcParams.update({'font.size': 14})

    ax1.plot(diameters, T_iabb_L2,
             color="black", marker="x", markevery=1, markersize=8, linestyle="", linewidth=2,
             label='IABB')

    ax1.plot(diameters, T_abb_L2,
             color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
             label='ABB')

    ax1.plot(diameters, y_1st,
             color="black", marker="", markevery=1, markersize=5, linestyle="--", linewidth=2,
             label=r'$\mathcal{O}(n)$ convergence')

    ax1.plot(diameters, y_2nd,
             color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
             label=r'$\mathcal{O}(n^2)$ convergence')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks(diameters)

    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # plt.title(f'Pipe within pipe Benchmark - Grid Convergence Study\n '
    #           r'$k$=' + f'{k} \t')
    plt.xlabel(r'lattice size [lu]', fontsize=18)
    plt.ylabel(r'$u_{x}: \; L_2 \, error \, norm $', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=1E-16)
    plt.legend()
    plt.grid(True,  which="both")
    # plt.grid(True)
    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()
    plt.close(fig)  # close the figure


###################################################################################################################
def see_ux_plot(d):
    print("------------------------------------ ux  PLOT ------------------------------------")
    fig_name = f'Conduction_between_plates_anal_vs_lbm_k{conductivity}_q{q}_effdiam_{effdiam}.png'

    # -------------------- make dummy plot --------------------
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(14, 8))

    axes = plt.gca()

    plt.plot(T_anal, y_anal + expected_wall_location,
             color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
             label=r'$analytical \, solution$')

    # plt.plot(u_anal, y_anal + len(y_grid) / 2,
    #          color="black", marker="o", markevery=1, markersize=5, linestyle=":", linewidth=2,
    #          label=r'$analytical \, solution$')

    # plt.plot(u_fd, y_fd + len(y_grid) / 2,
    #          color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
    #          label=r'$FD \, solution$')

    plt.plot(T_iabb_num_slice, y_grid,
             color="black", marker="x", markevery=1, markersize=7, linestyle="", linewidth=2,
             label=r'$LBM \, IBB$')

    plt.plot(T_abb_num_slice, y_grid,
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


    plt.title(f'IABB conductucion_between_plates:\n' +
              r'$ \nu = $' + f'{conductivity} ' + r'$D_{eff}=$' + f'{effdiam}, q={q}')

    plt.xlabel(r'$u_x$')
    plt.ylabel(r'$y$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    # plt.close(fig)  # close the figure
