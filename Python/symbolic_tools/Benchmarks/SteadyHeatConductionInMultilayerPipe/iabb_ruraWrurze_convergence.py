
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.ticker
from Benchmarks.PoiseuilleFlow.PoiseuilleAnal import calc_gx_in_pipe, OnePhasePoiseuilleAnalInPipe
# from Benchmarks.PoiseuilleFlow.pipe_poiseuille_plots import cntr_plot, slice_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeDirichlet

from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, get_r_from_xy
from DataIO.helpers import calc_L2, calc_mse
import pwd
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.VTIFile import VTIFile

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IABB_ruraWrurze')


eff_pipe_diam = np.array([30, 46, 60, 92, 120])
eff_cyl_diam = np.array([15, 23, 30, 46, 60])
conductivities = np.array([0.1])
kin_visc = 0.1

x0 = 64.  # center of the cylinder/pipe
y0 = 64.  # center of the cylinder/pipe

output_format = "png"

# prepare storage
n_diam = len(eff_pipe_diam)
n_conductivities = len(conductivities)

T_iabb_mse = np.zeros([n_conductivities, n_diam])
T_iabb_L2 = np.zeros([n_conductivities, n_diam])
T_abb_mse = np.zeros([n_conductivities, n_diam])
T_abb_L2 = np.zeros([n_conductivities, n_diam])


def read_data_from_LBM(case_folder):
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    T_num = vti_reader.get("T")
    [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    T_num_slice = T_num[:, :, 1]
    # ny, nx, nz = T_num.shape
    # uz_num_slice = uz_num[:, :, 1]
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return T_num_slice


def calculate_error_norms(_case_folder):
    T_num_slice = read_data_from_LBM(_case_folder)
    ny, nx = T_num_slice.shape

    r0 = eff_cyl_diam[d] / 2.  # inner radius
    r2 = eff_pipe_diam[d] / 2.  # outer radius
    r1 = (r0 + r2) / 2.  # interface between layers

    pwp = PipeWithinPipeDirichlet(r0, r1, r2, conductivities[k], conductivities[k], T0=11, T2=10)

    x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
    y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
    xx, yy = np.meshgrid(x_grid, y_grid)

    T_anal = np.zeros((ny, nx))
    r_anal = np.zeros((ny, nx))

    cuttoff_r2 = eff_pipe_diam[d] / 2. - 1
    cuttoff_r0 = eff_cyl_diam[d] / 2. + 1

    for i in range(ny):
        for j in range(nx):
            r = get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
            r_anal[i, j] = r
            T_anal[i, j] = pwp.get_temperature_r(r)
            if r < cuttoff_r0 or r > cuttoff_r2:
                T_anal[i, j] = np.nan

    not_nan_mask = ~np.isnan(T_anal)
    T_anal_masked = T_anal[not_nan_mask]
    T_num_slice_masked = T_num_slice[not_nan_mask]
    r_anal_masked = r_anal[not_nan_mask]

    cntr_plot(T_anal, T_num_slice, xx, yy, conductivities[k], eff_pipe_diam[d])

    # u_anal_masked = np.array([max(u_anal_masked)])
    # uz_num_slice_masked = np.array([max(uz_num_slice_masked)])
    mse = calc_mse(T_anal_masked, T_num_slice_masked)
    L2 = calc_L2(T_anal_masked, T_num_slice_masked)
    return mse, L2


for k in range(n_conductivities):
    for d in range(n_diam):
        case_folder = f'iabb_ruraWrurze_Dirichlet_Cumulants_k_{conductivities[k]}_nu_{kin_visc}_effdiam_{eff_pipe_diam[d]}'

        T_iabb_mse[k, d], T_iabb_L2[k, d] = calculate_error_norms(os.path.join(main_folder, case_folder))

        print(f"uz_ibb_mse={T_iabb_mse[k, d]:.2e} for k{conductivities[k]}_effdiam_{eff_pipe_diam[d]}")
        print(f"uz_ibb_L2={T_iabb_L2[k, d]:.2e} for v{conductivities[k]}_effdiam_{eff_pipe_diam[d]}")

        case_folder = f'abb_ruraWrurze_Dirichlet_Cumulants_k_{conductivities[k]}_nu_{kin_visc}_effdiam_{eff_pipe_diam[d]}'
        T_abb_mse[k, d], T_abb_L2[k, d] = calculate_error_norms(os.path.join(main_folder, case_folder))

        print(f"uz_bb_mse={T_abb_mse[k, d]:.2e} for k{conductivities[k]}_effdiam_{eff_pipe_diam[d]}")
        print(f"uz_bb_L2={T_abb_L2[k, d]:.2e} for k{conductivities[k]}_effdiam_{eff_pipe_diam[d]}")


def make_plot_for_given_conductivity(_k):
    print("------------------------------------ Convergence  PLOT ------------------------------------")
    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/grid_convergence_iabb_ruraWrurze_anal_vs_lbm_k{eat_dots_for_texmaker(conductivities[_k])}.{output_format}'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    initial_error_05st = 0.070
    y_05st = np.sqrt(eff_pipe_diam.min())*initial_error_05st/np.sqrt(eff_pipe_diam)

    # initial_error_1st = 0.18
    initial_error_1st = 1.15*max(np.concatenate((T_iabb_L2[_k, :], T_abb_L2[_k, :])))
    y_1st = eff_pipe_diam.min()*initial_error_1st/eff_pipe_diam
    # initial_error_2nd = 0.05
    initial_error_2nd = 0.85 * min((max(T_iabb_L2[_k, :]), max(T_abb_L2[_k, :])))
    # initial_error_2nd = 0.85 * max((max(uz_ibb_L2[_v, :]), max(uz_bb_L2[_v, :])))
    y_2nd = eff_pipe_diam.min()*eff_pipe_diam.min()*initial_error_2nd/(eff_pipe_diam*eff_pipe_diam)

    fig1, ax1 = plt.subplots(figsize=(14, 8))
    plt.rcParams.update({'font.size': 14})

    ax1.plot(eff_pipe_diam, T_iabb_L2[_k, :],
             color="black", marker="x", markevery=1, markersize=8, linestyle="", linewidth=2,
             label='IABB')

    ax1.plot(eff_pipe_diam, T_abb_L2[_k, :],
             color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
             label='ABB')

    ax1.plot(eff_pipe_diam, y_1st,
             color="black", marker="", markevery=1, markersize=5, linestyle="--", linewidth=2,
             label=r'$\mathcal{O}(n)$ convergence')

    ax1.plot(eff_pipe_diam, y_2nd,
             color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
             label=r'$\mathcal{O}(n^2)$ convergence')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks(eff_pipe_diam)

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


def make_plot_for_all_viscosities():
    print("------------------------------------ Convergence  PLOT 2 ------------------------------------")

    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/grid_convergence_PipePoiseuille_all.{output_format}'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    initial_error_05st = 0.070
    y_05st = np.sqrt(eff_pipe_diam.min())*initial_error_05st/np.sqrt(eff_pipe_diam)

    # initial_error_1st = 0.18
    # initial_error_1st = 1.15*max(np.concatenate((uz_ibb_L2[:, :], uz_bb_L2[:, :])))
    initial_error_1st = 1.15 * max(np.concatenate([np.concatenate(T_iabb_L2.tolist()), np.concatenate(T_abb_L2.tolist())]))
    y_1st = eff_pipe_diam.min()*initial_error_1st/eff_pipe_diam
    # initial_error_2nd = 0.05
    initial_error_2nd = 0.85*min((max(np.concatenate(T_iabb_L2.tolist())), max(np.concatenate(T_abb_L2.tolist()))))
    y_2nd = eff_pipe_diam.min()*eff_pipe_diam.min()*initial_error_2nd/(eff_pipe_diam*eff_pipe_diam)



    fig1, ax1 = plt.subplots(figsize=(14, 8))
    plt.rcParams.update({'font.size': 14})

    markers = ["x", ">", "o", "v", "d"]
    for k_idx in range(len(conductivities)):
        ax1.plot(eff_pipe_diam, T_iabb_L2[k_idx, :],
                 color="black", marker=markers[k_idx], markevery=1, markersize=8, linestyle="", linewidth=2,
                 label=r'IBB: $\nu$ =' + f'{conductivities[k_idx]}')

    # ax1.plot(effdiams, uz_bb_L2[k_visc_ind + 1, :],
    #          color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
    #          label=r'IBB - $\nu$' + f'{kin_viscs[k_visc_ind + 1 ]}')

    ax1.plot(eff_pipe_diam, y_1st,
             color="black", marker="", markevery=1, markersize=5, linestyle="--", linewidth=2,
             label=r'$\mathcal{O}(n)$ convergence')

    ax1.plot(eff_pipe_diam, y_2nd,
             color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
             label=r'$\mathcal{O}(n^2)$ convergence')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks(eff_pipe_diam)

    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # plt.title(f'Pipe within pipe Benchmark - Grid Convergence Study\n '
    #           r'$k$=' + f'{k} \t')
    plt.xlabel(r'lattice size [lu]', fontsize=20)
    plt.ylabel(r'$u_{x}: \; L_2 \, error \, norm $', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='minor', labelsize=1E-16)
    plt.legend()
    plt.grid(True,  which="both")
    # plt.grid(True)
    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()
    plt.close(fig)  # close the figure


for k in range(n_conductivities):
    make_plot_for_given_conductivity(k)

# make_plot_for_all_viscosities()
