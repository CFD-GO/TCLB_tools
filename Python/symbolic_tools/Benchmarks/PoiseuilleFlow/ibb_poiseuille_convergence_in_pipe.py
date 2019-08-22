
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.ticker
from Benchmarks.PoiseuilleFlow.PoiseuilleAnal import calc_gx_in_pipe, OnePhasePoiseuilleAnalInPipe
from Benchmarks.PoiseuilleFlow.pipe_poiseuille_plots import cntr_plot, slice_plot
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
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IBB_PoiseuillePipe')


rho = 1.
kin_viscs = np.array([0.1, 0.1666666, 0.01, 0.001])
effdiams = np.array([15, 23, 31, 45, 61, 93, 121])

x0 = 63.5  # center of the pipe
y0 = 63.5

uc = 0.01

output_format = "pdf"

# prepare storage
n_diam = len(effdiams)
n_kin_visc = len(kin_viscs)
uz_ibb_mse = np.zeros([n_kin_visc, n_diam])
uz_ibb_L2 = np.zeros([n_kin_visc, n_diam])
uz_bb_mse = np.zeros([n_kin_visc, n_diam])
uz_bb_L2 = np.zeros([n_kin_visc, n_diam])


def read_data_from_LBM(case_folder):
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    T_num = vti_reader.get("T")
    [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    ny, nx, nz = T_num.shape

    uz_num_slice = uz_num[:, :, 1]
    # T_num_slice = T_num[:, :, 1]
    #
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return uz_num_slice


def calculate_error_norms(_case_folder):
    uz_num_slice = read_data_from_LBM(_case_folder)
    ny, nx = uz_num_slice.shape

    gx = calc_gx_in_pipe(uc, rho, kin_viscs[k], D=effdiams[d])
    poiseuilleAnal = OnePhasePoiseuilleAnalInPipe(gx=gx, nu=kin_viscs[k], D=effdiams[d])

    x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
    y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
    xx, yy = np.meshgrid(x_grid, y_grid)

    u_anal = np.zeros((ny, nx))
    r_anal = np.zeros((ny, nx))
    r_cutoff = effdiams[d] / 2. - 1

    for i in range(ny):
        for j in range(nx):
            r = get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
            r_anal[i, j] = r
            u_anal[i, j] = poiseuilleAnal.get_u_profile(r)
            if r > r_cutoff:
                u_anal[i, j] = np.nan

    not_nan_mask = ~np.isnan(u_anal)
    u_anal_masked = u_anal[not_nan_mask]
    uz_num_slice_masked = uz_num_slice[not_nan_mask]
    r_anal_masked = r_anal[not_nan_mask]

    # cntr_plot(u_anal, uz_num_slice, kin_viscs[k], effdiams[d], nx, xx, yy)
    #
    # u_anal_masked = np.array([max(u_anal_masked)])
    # uz_num_slice_masked = np.array([max(uz_num_slice_masked)])
    mse = calc_mse(u_anal_masked, uz_num_slice_masked)
    L2 = calc_L2(u_anal_masked, uz_num_slice_masked)
    return mse, L2


for k in range(n_kin_visc):
    for d in range(n_diam):
        case_folder = f'ibb_PoiseuillePipe_CM_HIGHER_nu_{kin_viscs[k]}_effdiam_{effdiams[d]}'
        uz_ibb_mse[k, d], uz_ibb_L2[k, d] = calculate_error_norms(os.path.join(main_folder, case_folder))

        print(f"uz_ibb_mse={uz_ibb_mse[k,d]:.2e} for v{kin_viscs[k]}_effdiam_{effdiams[d]}")
        print(f"uz_ibb_L2={uz_ibb_L2[k,d]:.2e} for v{kin_viscs[k]}_effdiam_{effdiams[d]}")

        case_folder = f'bb_PoiseuillePipe_CM_HIGHER_nu_{kin_viscs[k]}_effdiam_{effdiams[d]}'
        uz_bb_mse[k, d], uz_bb_L2[k, d] = calculate_error_norms(os.path.join(main_folder, case_folder))

        print(f"uz_bb_mse={uz_bb_mse[k,d]:.2e} for v{kin_viscs[k]}_effdiam_{effdiams[d]}")
        print(f"uz_bb_L2={uz_bb_L2[k,d]:.2e} for v{kin_viscs[k]}_effdiam_{effdiams[d]}")


def make_plot_for_given_viscosity(_v):
    print("------------------------------------ Convergence  PLOT ------------------------------------")
    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/grid_convergence_PipePoiseuille_anal_vs_lbm_nu{eat_dots_for_texmaker(kin_viscs[_v])}.{output_format}'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    initial_error_05st = 0.070
    y_05st = np.sqrt(effdiams.min())*initial_error_05st/np.sqrt(effdiams)

    # initial_error_1st = 0.18
    initial_error_1st = 1.15*max(np.concatenate((uz_ibb_L2[_v, :], uz_bb_L2[_v, :])))
    y_1st = effdiams.min()*initial_error_1st/effdiams
    # initial_error_2nd = 0.05
    initial_error_2nd = 0.85 * min((max(uz_ibb_L2[_v, :]), max(uz_bb_L2[_v, :])))
    # initial_error_2nd = 0.85 * max((max(uz_ibb_L2[_v, :]), max(uz_bb_L2[_v, :])))
    y_2nd = effdiams.min()*effdiams.min()*initial_error_2nd/(effdiams*effdiams)

    fig1, ax1 = plt.subplots(figsize=(14, 8))
    plt.rcParams.update({'font.size': 14})

    ax1.plot(effdiams, uz_ibb_L2[_v, :],
             color="black", marker="x", markevery=1, markersize=8, linestyle="", linewidth=2,
             label='IBB')

    ax1.plot(effdiams, uz_bb_L2[_v, :],
             color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
             label='BB')

    ax1.plot(effdiams, y_1st,
             color="black", marker="", markevery=1, markersize=5, linestyle="--", linewidth=2,
             label=r'$\mathcal{O}(n)$ convergence')

    ax1.plot(effdiams, y_2nd,
             color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
             label=r'$\mathcal{O}(n^2)$ convergence')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks(effdiams)

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
    y_05st = np.sqrt(effdiams.min())*initial_error_05st/np.sqrt(effdiams)

    # initial_error_1st = 0.18
    # initial_error_1st = 1.15*max(np.concatenate((uz_ibb_L2[:, :], uz_bb_L2[:, :])))
    initial_error_1st = 1.15 * max(np.concatenate([np.concatenate(uz_ibb_L2.tolist()), np.concatenate(uz_bb_L2.tolist())]))
    y_1st = effdiams.min()*initial_error_1st/effdiams
    # initial_error_2nd = 0.05
    initial_error_2nd = 0.85*min((max(np.concatenate(uz_ibb_L2.tolist())),max(np.concatenate(uz_bb_L2.tolist()))))
    y_2nd = effdiams.min()*effdiams.min()*initial_error_2nd/(effdiams*effdiams)



    fig1, ax1 = plt.subplots(figsize=(14, 8))
    plt.rcParams.update({'font.size': 14})

    markers = ["x", ">", "o", "v", "d"]
    for k_visc_ind in range(len(kin_viscs)):
        ax1.plot(effdiams, uz_ibb_L2[k_visc_ind, :],
                 color="black", marker=markers[k_visc_ind], markevery=1, markersize=8, linestyle="", linewidth=2,
                 label=r'IBB: $\nu$ =' + f'{kin_viscs[k_visc_ind]}')

    # ax1.plot(effdiams, uz_bb_L2[k_visc_ind + 1, :],
    #          color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
    #          label=r'IBB - $\nu$' + f'{kin_viscs[k_visc_ind + 1 ]}')

    ax1.plot(effdiams, y_1st,
             color="black", marker="", markevery=1, markersize=5, linestyle="--", linewidth=2,
             label=r'$\mathcal{O}(n)$ convergence')

    ax1.plot(effdiams, y_2nd,
             color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
             label=r'$\mathcal{O}(n^2)$ convergence')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks(effdiams)

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


for v in range(n_kin_visc):
    make_plot_for_given_viscosity(v)

make_plot_for_all_viscosities()
