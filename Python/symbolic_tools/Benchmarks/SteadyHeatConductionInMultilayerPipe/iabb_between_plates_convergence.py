
import numpy as np
import os
from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import HeatConductionBetweenTwoPlates
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, calc_mse, calc_L2, strip_folder_name, eat_dots_for_texmaker
import pwd
from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot, slice_plot
from DataIO.VTIFile import VTIFile

from matplotlib.ticker import LinearLocator, FormatStrFormatter
from DataIO.VTIFile import VTIFile
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.pylab as pylab
rho = 1
conductivity = '0.1'
mu = rho * float(conductivity)

diameters = np.array([15, 31, 63, 127])
qs = [0.25, 0.5, 0.75]

collision_type = 'CM_HIGHER'

# -------- numerical solution ---------------


wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IABB_plates')


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

        iabb_case_folder = f"iabb_plates_Dirichlet_{collision_type}_k_{conductivity}_effdiam_{effdiam}"
        T_iabb_num_slice = read_ux(iabb_case_folder)
        abb_case_folder = f"abb_plates_Dirichlet_{collision_type}_k_{conductivity}_effdiam_{effdiam}"
        T_abb_num_slice = read_ux(abb_case_folder)

        # T_iabb_num_slice = delete_unphysical_data_from_wall_nodes(T_iabb_num_slice)
        # T_abb_num_slice = delete_unphysical_data_from_wall_nodes(T_abb_num_slice)
        # y_grid = delete_unphysical_data_from_wall_nodes(y_grid)


        hcbp = HeatConductionBetweenTwoPlates(T0=11, T2=10, Heff=effdiam, y0=expected_wall_location)

        ny, nx = T_iabb_num_slice.shape
        x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
        y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
        xx, yy = np.meshgrid(x_grid, y_grid)

        T_anal = np.zeros((ny, nx))

        cuttoff_y2 = effdiam - 2
        cuttoff_y0 = 0 + 2
        for i in range(ny):
            for j in range(nx):
                y = yy[i][j]
                T_anal[i, j] = hcbp.get_T_profile(y)
                if i < cuttoff_y0 or i > cuttoff_y2:
                    T_anal[i, j] = np.nan

        not_nan_mask = ~np.isnan(T_anal)
        T_anal_masked = T_anal[not_nan_mask]
        T_iabb_num_slice_masked = T_iabb_num_slice[not_nan_mask]
        T_abb_num_slice_masked = T_abb_num_slice[not_nan_mask]
        yy_masked = yy[not_nan_mask]

        T_iabb_mse[d] = calc_mse(T_anal_masked, T_iabb_num_slice_masked)
        T_iabb_L2[d] = calc_L2(T_anal_masked, T_iabb_num_slice_masked)
        T_abb_mse[d] = calc_mse(T_anal_masked, T_abb_num_slice_masked)
        T_abb_L2[d] = calc_L2(T_anal_masked, T_abb_num_slice_masked)


        # cntr_plot(T_anal, T_num_slice, xx, yy, conductivity, effdiam, title=iabb_case_folder)

        # # 2D clipy_grid
        # yy_masked = yy_masked[:,  int(nx / 2)]
        T_anal_slice = T_anal[:, int(nx / 2)]
        T_iabb_num_slice = T_iabb_num_slice[:, int(nx / 2)]
        not_nan_mask = ~np.isnan(T_anal_slice)
        T_anal_masked = T_anal_slice[not_nan_mask]
        T_iabb_num_slice_masked = T_iabb_num_slice[not_nan_mask]
        y_masked = y_grid[not_nan_mask]

        slice_plot(T_anal_masked, T_iabb_num_slice_masked, y_masked, title=f'{iabb_case_folder}_q{q}')

        print(f"T_iabb_mse={T_iabb_mse[d]:.2e} T_abb_mse={T_abb_mse[d]:.2e} \t for k{conductivity}_q{q}_effdiam_{effdiam}")
        print(f"T_iabb_L2={T_iabb_L2[d]:.2e} T_abb_L2={T_abb_L2[d]:.2e} \tfor k{conductivity}_q{q}_effdiam_{effdiam}")

    print("------------------------------------ Convergence  PLOT ------------------------------------")
    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/grid_convergence_conduction_between_plates_k{eat_dots_for_texmaker(conductivity)}_q{eat_dots_for_texmaker(q)}.png'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    initial_error_05st = 0.070
    y_05st = np.sqrt(diameters.min())*initial_error_05st/np.sqrt(diameters)


    initial_error_1st = 1.15*max(np.concatenate((T_iabb_L2, T_abb_L2)))
    # initial_error_1st = 1.15*max((T_iabb_L2))
    # initial_error_1st = 1.15 * max((T_abb_L2))
    y_1st = diameters.min()*initial_error_1st/diameters

    initial_error_2nd = 0.85*min((max(T_iabb_L2), max(T_abb_L2)))
    # initial_error_2nd = 0.85 * (max(T_iabb_L2))
    # initial_error_2nd = 0.85 * (max(T_abb_L2))
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
    plt.ylabel(r'$T: \; L_2 \, error \, norm $', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=1E-16)
    plt.legend()
    plt.grid(True,  which="both")
    # plt.grid(True)
    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()
    plt.close(fig)  # close the figure
