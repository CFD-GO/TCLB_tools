import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, delete_unphysical_data_from_wall_nodes
import pwd
from Benchmarks.TransientHeatConductionIn1DCompositeSlab.TransientCHTin1DCompositeSlab_easy import Solver
from DataIO.VTIFile import VTIFile

import matplotlib.pylab as pylab


# Notation
# Follows Fig.1. from 'Phase interface effects in the total enthalpy-based lattice Boltzmann model for solid–liquid phase change'
T_left = 1
T_right = 0
k_right = 0.1
cp_right = 1
# R_cp = 1./4  # R_cp = cp_left/cp_right
R_k = 8  # R_k = k_left/k_right
R_cps = [1./16, 1./4, 1, 4, 16]  # R_cp = cp_left/cp_right
# R_ks = [1./8, 1, 8]  # R_k = k_left/k_right


# gamma = 1  # stability enhancement
str_iteration = "00010001"
nodes_to_strip_per_side = 925  # number of extra BC nodes in domain per side

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_cht_lem_MRT')


def read_data_from_LBM(case_folder, iteration_of_interest):
    # iteration_of_interest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, iteration_of_interest, extension='.vti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=False)
    fs_num = vti_reader.get("FractionSuspected")
    # [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    # ny, nx, nz = T_num.shape
    # uz_num_slice = uz_num[:, :, 1]

    fs_num_slice = fs_num[int(fs_num.shape[0]/2), :]
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return fs_num_slice

def prepare_LBM_and_anal_data_pair(R_k, R_cp, T_left, T_right, k_right, cp_right):
    k_left = R_k*k_right
    cp_left = R_cp*cp_right
    solver = Solver(T_left, T_right, k_left, k_right, cp_left, cp_right)

    folder = f"lem_R_k_{R_k:.2e}_R_cp_{R_cp:.2e}_k_right_{k_right}_cp_right_{cp_right}_size_2048lu"
    case_folder = os.path.join(main_folder, folder)

    temperature_lbm = read_data_from_LBM(case_folder, str_iteration)
    # The Dirichlet BC is imposed as wet node BC (equilibrium) in LBM simulation
    temperature_lbm = delete_unphysical_data_from_wall_nodes(temperature_lbm, nodes_to_strip_per_side)
    effective_domain = temperature_lbm.shape[0]

    # -------- analytical solution ---------------
    dx = (1+1) / effective_domain  # (physical domain is from -1 to 1)/(LBM domain length)
    x_grid = np.linspace(start=-1+dx/2., stop=1-dx/2., num=effective_domain, endpoint=True)

    k_right_lbm = k_right  # for simplicity, diffusivity_physical matches diffusivity_LBM
    time_step = (dx*dx)*k_right_lbm/k_right  # describes how many second is in one lbm timestep
    time_spot = int(str_iteration) * time_step
    temperature_anal = [solver.calc_transient_T_profile(time_spot, x_i) for x_i in x_grid]

    return x_grid, temperature_lbm, temperature_anal, time_spot

def prepare_plot():
    if not os.path.exists('plots'):
        os.makedirs('plots')

    plt.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(14, 8))
    axes = fig.gca()

    # ------ format y axis ------ #
    # yll = y.min()
    # yhl = y.max()
    # axes.set_ylim([yll, yhl])

    # axes.set_yticks(np.linspace(yll, yhl, 5))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    axes.set_ylim([0, 1])
    axes.set_yticks([0, 0.25, 0.5, 0.75, 1])
    # axes.yaxis.set_major_formatter(xfmt)

    # plt.yscale('log')


    # ------ format x axis ------ #
    axes.set_xlim([-1.1, 1.1])
    # plt.xlim(x1-0.5, x2+0.5)
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    axes.set_title(
        f'T across interface\n '
        r"$R_{k}=$" + f"{R_k}")

    axes.set_xlabel(r'$x$')
    axes.set_ylabel(r'$Temperature$')

    return axes, fig


def add_data_to_plot(axes, marker, x_grid, temperature_lbm, temperature_anal, time_spot, R_cp):
    axes.plot(x_grid, temperature_anal,
             color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
             label=f'analytical R_cp {R_cp} - t={time_spot:.2f}[s]')

    axes.plot(x_grid, temperature_lbm,
             color="black", marker=marker, markevery=7, markersize=10, linestyle=":", linewidth=2,
             label=f'LBM R_cp {R_cp} - at t={time_spot:.2f}[s]')

    axes.legend()
    axes.grid()


axes, fig = prepare_plot()

plt_markers = ["<", "v", ">", "x", "o"]
for R_cp, plt_marker in zip(R_cps, plt_markers):
    x_grid, temperature_lbm, temperature_anal, time_spot = \
        prepare_LBM_and_anal_data_pair(R_k, R_cp, T_left, T_right, k_right, cp_right)
    add_data_to_plot(axes, plt_marker, x_grid, temperature_lbm, temperature_anal, time_spot, R_cp)

R_k_tex = eat_dots_for_texmaker(R_k)
fig_name = f'plots/TransientHeatConductionInCompositeSlab_Rk{R_k_tex}.png'
fig.savefig(fig_name, bbox_inches='tight')

