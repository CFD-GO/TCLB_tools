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
R_cp = 1  # R_cp = cp_left/cp_right
R_k = 1.  # R_k = k_left/k_right
# R_cp = [1./16, 1./4, 1, 4, 16]  # R_cp = cp_left/cp_right
# R_k = [1./8, 1, 8]  # R_k = k_left/k_right
k_left = R_k*k_right
cp_left = R_cp*cp_right

gamma = 1  # stability enhancement

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

    fs_num_slice = fs_num[int(fs_num.shape[1]/2), :]
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return fs_num_slice


# folder_hard = f"lem_R_k_1.00e+00_R_cp_1.00e+00_k_right_0.1_cp_right_1_size_64lu"
folder = f"lem_R_k_{R_k:.2e}_R_cp_{R_cp:.2e}_k_right_{k_right}_cp_right_{cp_right}_size_64lu"
case_folder = os.path.join(main_folder, folder)

# str_iteration = "00000501"
str_iteration = "00001001"
temperature_lbm = read_data_from_LBM(case_folder, str_iteration)

# The Dirichlet BC is imposed as wet node BC (equilibrium) in LBM simulation
nodes_to_strip_per_side = 1  # number of extra BC nodes in domain per side
temperature_lbm = delete_unphysical_data_from_wall_nodes(temperature_lbm, nodes_to_strip_per_side)

xSIZE = temperature_lbm.shape[0]
effective_domain = xSIZE

# -------- analytical solution ---------------
solver = Solver(T_left, T_right, k_left, k_right, cp_left, cp_right)
dx = (1+1) / effective_domain  # (physical domain is from -1 to 1)/(LBM domain length)
x_grid = np.linspace(start=-1+dx/2., stop=1-dx/2., num=xSIZE, endpoint=True) # + 0.5


k_right_lbm = k_right
time_step = (dx*dx)*k_right_lbm/k_right  # describes how many second is in one lbm timestep

#  time_spot = 1200*1/(60*60) = 0.33333333
time_spot = (int(str_iteration)-0) * time_step
# time_spot += 1.615
temperature_anal = [solver.calc_transient_T_profile(time_spot, x_i) for x_i in x_grid]


# -------------------- make dummy plot --------------------
if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/TransientHeatConductionInCompositeSlab.png'

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

plt.plot(x_grid, temperature_anal,
         color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label=f'transient analytical solution - at time={time_spot}[s]')

plt.plot(x_grid, temperature_lbm,
         color="black", marker="v", markevery=7, markersize=10, linestyle=":", linewidth=2,
         label=f'transient LBM solution - at time={time_spot}[s]')


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


plt.title(f'T across interface\n '
          r"$R_{cp}=$"+ f"{R_cp}" + " ; " + r"$R_{k}=$" + f"{R_k}"
          # r'$x_{1}$=' + f'{x1}' + '\t' + r'$x_{2}$=' + f'{x2}'
          # f'; \t'
          # r'$x_{step}$' + f'={step:.4f}'
          )
plt.xlabel(r'$x$')
plt.ylabel(r'$Temperature$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()
