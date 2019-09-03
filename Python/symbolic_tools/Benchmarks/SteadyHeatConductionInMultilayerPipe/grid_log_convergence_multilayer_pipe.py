

from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeDirichlet
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from DataIO.VTIFile import VTIFile
import os
import pwd
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

basic_lattice_size = 32
gauges = np.array([1, 2, 4, 8])
lattices = gauges * basic_lattice_size

k_inner = 0.1  # k_inner is fixed
k_outer = 0.01  # 0.01, 0.001

def get_t_mse(folder):
    n = len(gauges)
    T_mse = np.zeros(n)
    T_L2 = np.zeros(n)

    for g in range(n):
        filepath_vtk = os.path.join(folder,
                                    f'k_outer_{k_outer}_size_{lattices[g]}lu',
                                    f'k_outer_{k_outer}_size_{lattices[g]}lu_VTK_P00_00500000.vti')
        vti_reader = VTIFile(filepath_vtk)

        T_num = vti_reader.get("T")
        ySIZE, xSIZE = T_num.shape
        assert ySIZE == xSIZE == lattices[g]
        assert xSIZE % 1 == 0

        r0 = gauges[g] * (8 / 2)  # inner radius
        r2 = gauges[g] * (30 / 2)  # outer radius

        abb_correction = 0.5
        if 'abb_scheme' in filepath_vtk:
            r0 -= abb_correction
            r2 += abb_correction

        r1 = gauges[g] * (20 / 2)  # interface between layers
        x0 = lattices[g] / 2  # center of the pipe
        y0 = lattices[g] / 2

        # ----------------------- compute anal solution ---------------------------

        pwp = PipeWithinPipeDirichlet(r0, r1, r2, k_inner, k_outer, T0=0, T2=1)

        x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
        y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
        xx, yy = np.meshgrid(x_grid, y_grid)
        T_anal = np.zeros((ySIZE, xSIZE))

        for i in range(ySIZE):
            # print(f"=== Doing i/ny: {i}/{ny}  ===")
            for j in range(xSIZE):
                # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
                r = pwp.get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
                T_anal[i][j] = pwp.get_temperature_r(r)

                # if r < r0 or r > r2:
                #     T_anal[i][j] = 0
                #     T_num[i][j] = 0

        # 2d clip
        T_anal = T_anal[:, int(xSIZE / 2)]
        T_num = T_num[:, int(xSIZE / 2)]

        T_L2[g] = np.sqrt(
            np.sum((T_anal - T_num) * (T_anal - T_num))
            / np.sum(T_anal * T_anal))  # Eq. 4.57


        T_mse[g] = np.sum((T_anal - T_num) * (T_anal - T_num)) / len(T_anal)
        print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
        print(f"T_L2={T_L2[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
    return T_L2
    # return T_mse


home = pwd.getpwuid(os.getuid()).pw_dir
T_err_EQ = get_t_mse(os.path.join(home, 'DATA_FOR_PLOTS', 'batch_ruraWrurze_variable_k', 'eq_scheme'))
T_err_ABB = get_t_mse(os.path.join(home, 'DATA_FOR_PLOTS', 'batch_ruraWrurze_variable_k', 'abb_scheme'))


print("------------------------------------ PLOT ------------------------------------")
fig_name = f'pipe_within_pipe_grid_convergence_k_inner{k_inner}_k_outer{k_outer}.png'

initial_error_05st = 0.070
y_05st = np.sqrt(lattices.min())*initial_error_05st/np.sqrt(lattices)
initial_error_1st = 0.2
# initial_error_1st = 0.044
# initial_error_1st = 4.7e-4  # mse
y_1st = lattices.min()*initial_error_1st/lattices
initial_error_2nd = 0.01
# initial_error_2nd = 1.0e-5  # mse
y_2nd = lattices.min()*lattices.min()*initial_error_2nd/(lattices*lattices)


fig1, ax1 = plt.subplots(figsize=(14, 8))
plt.rcParams.update({'font.size': 14})

ax1.plot(lattices, T_err_EQ,
         color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
         label='Equilibrium scheme')

ax1.plot(lattices, T_err_ABB,
         color="black", marker=">", markevery=1, markersize=5, linestyle="", linewidth=2,
         label='Anti-Bounce-Back scheme')

# ax1.plot(lattices, y_05st,
#          color="black", marker="", markevery=1, markersize=5, linestyle=":", linewidth=2,
#          label=r'$\mathcal{O}(\sqrt{n})$ convergence')

ax1.plot(lattices, y_1st,
         color="black", marker="", markevery=1, markersize=5, linestyle="--", linewidth=2,
         label=r'$\mathcal{O}(n)$ convergence')

ax1.plot(lattices, y_2nd,
         color="black", marker="", markevery=1, markersize=5, linestyle="-", linewidth=2,
         label=r'$\mathcal{O}(n^2)$ convergence')


ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xticks(lattices)


ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.title(f'Pipe within pipe Benchmark - Grid Convergence Study\n '
          r'$k_{inner}$=' + f'{k_inner} \t;\t' 
          r'$k_{outer}$=' + f'{k_outer} '
          # r'$x_{range}$' + f'={range}'
          # f'; \t'
          # r'$x_{step}$' + f'={step:.4f}'
          )
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
# plt.close(fig)  # close the figure

print("bye")



