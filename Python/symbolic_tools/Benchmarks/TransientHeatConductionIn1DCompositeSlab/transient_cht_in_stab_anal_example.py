import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker
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
R_cp = 1.  # R_cp = cp_left/cp_right
R_k = 1.  # R_k = k_left/k_right
# R_cp = [1./16, 1./4, 1, 4, 16]  # R_cp = cp_left/cp_right
# R_k = [1./8, 1, 8]  # R_k = k_left/k_right
k_left = R_k * k_right
cp_left = R_cp * cp_right

# -------- numerical solution ---------------


# -------- analytical solution ---------------
solver = Solver(T_left, T_right, k_left, k_right, cp_left, cp_right)

x_grid = np.linspace(-1, 1, num=1000, endpoint=True)
time_spot = 1
ts = [solver.calc_transient_T_profile(time_spot, x_i) for x_i in x_grid]

# -------------------- make dummy plot --------------------
if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/TransientHeatConductionInCompositeSlab.png'

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

plt.plot(x_grid, ts,
         color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label=f'transient analytical solution - at time={time_spot}')

# plt.plot(x_grid, FS_lbm,
#          color="black", marker="v", markevery=10, markersize=10, linestyle="--", linewidth=2,
#          label=f'transient lbm solution - at time={time_spot}')

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

# plt.xlim(x1-0.5, x2+0.5)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'T across interface\n '
          r"$R_{cp}=$" + f"{R_cp}" + " ; " + r"$R_{k}=$" + f"{R_k}"
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
