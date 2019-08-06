
from DataIO.VTIFile import VTIFile
import os
import pandas as pd
import pwd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import collections
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotBarman3D')


Re = 10
Pr = 10


# README:
# size 1 = LB mesh: 1000x150x3 -- sample 150 y points from toolbox
# size 2 = LB mesh: 2000x300x6 -- sample 300 y points from toolbox
# size 4 = LB mesh: 4000x600x12 -- sample 600 y points from toolbox


def read_data_from_LBM(case_folder):
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    T_num = vti_reader.get("T")
    [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    ny, nx, nz = T_num.shape

    ux_num_slice = ux_num[:, :, 1]
    T_num_slice = T_num[:, :, 1]

    y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return T_num_slice, ux_num_slice, y


def read_data_from_toolbox(mesh_type, ny):
    filepath_csv = os.path.join(main_folder, 'from_toolbox', f'Temperature_{ny}_{mesh_type}_Re_{Re}_Pr_{Pr}.csv')
    T_toolbox = np.loadtxt(filepath_csv, delimiter=',', skiprows=1)

    filepath_csv = os.path.join(main_folder, 'from_toolbox', f'ux_{ny}_{mesh_type}_Re_{Re}_Pr_{Pr}.csv')
    ux_toolbox = np.loadtxt(filepath_csv, delimiter=',', skiprows=1)
    y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return T_toolbox, ux_toolbox, y



print("data reading complete, lets plot it!")


cross_sections = collections.OrderedDict()
cross_sections['A'] = 270
cross_sections['B'] = 300
cross_sections['C'] = 330
cross_sections['D'] = 390
cross_sections['E'] = 450
cross_sections['F'] = 960
cross_sections['G'] = 999

# mesh='clmax_1'
mesh='nodes_787776'


T_qs, ux_qs, y_qs = read_data_from_toolbox('clmax_1', ny=int(150*4))
T_qs_struct, ux_qs_struct, y_qs_struct = read_data_from_toolbox('nodes_787776', ny=int(150*2))


# scaling_type = "keepU"
scaling_type = "keep_nu_and_k"
collision_kernel = 'CM_HIGHER'

T_lb_slice1, ux_lb1, y_lb1 = read_data_from_LBM(
    os.path.join(main_folder, f'batch_HotKarman3D_{collision_kernel}', f'{scaling_type}_sizer_{1}_Re{Re}_Pr{Pr}'))

# T_lb_slice2, ux_lb2, y_lb2 = read_data_from_LBM(
#     os.path.join(main_folder, f'{scaling_type}_sizer_{2}_Re{Re}_Pr{Pr}',
#                  f'HotKarman3D_template_sizer_{2}_Re{Re}_Pr{Pr}_VTK_P00_05600000.pvti'))  # Pr10

T_lb_slice2, ux_lb2, y_lb2 = read_data_from_LBM(
    os.path.join(main_folder, f'batch_HotKarman3D_{collision_kernel}', f'{scaling_type}_sizer_{2}_Re{Re}_Pr{Pr}'))  # Pr100 Pr1000


# T_lb_slice3, ux_lb3, y_lb3 = read_data_from_LBM(
#     os.path.join(main_folder, f'{scaling_type}_sizer_{4}_Re{Re}_Pr{Pr}',
#                  f'HotKarman3D_template_sizer_{4}_Re{Re}_Pr{Pr}_VTK_P00_00980000.pvti'))

T_lb_slice3, ux_lb3, y_lb3 = read_data_from_LBM(
    os.path.join(main_folder, f'batch_HotKarman3D_{collision_kernel}', f'{scaling_type}_sizer_{4}_Re{Re}_Pr{Pr}'))

for cross_section, x_cut in cross_sections.items():
    title = f'Temperature LBM vs FEM \n cross-section {cross_section}'
    title = ''  # skip title for latex report
    fig_name = f'plots/T_LBM_vs_ToolBox_{mesh}_Re{Re}_Pr{Pr}_cross_section_{cross_section}.pdf'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    axes = plt.gca()

    plt.plot(T_lb_slice1[:, x_cut*1], y_lb1,
             color="black", marker="", markevery=3, markersize=5, linestyle=":", linewidth=2,
             label='LBM - lattice I')

    plt.plot(T_lb_slice2[:, x_cut*2], y_lb2,
             color="black", marker="", markevery=4, markersize=5, linestyle="-.", linewidth=2,
             label='LBM - lattice II')

    plt.plot(T_lb_slice3[:, x_cut*4], y_lb3,
             color="black", marker="", markevery=5, markersize=5, linestyle="--", linewidth=2,
             label='LBM - lattice III')

    idx = list(cross_sections.keys()).index(cross_section)
    # plt.plot(T_qs[:, idx], y_qs,
    #          color="black", marker=">", markevery=5, markersize=5, linestyle="-", linewidth=2,
    #          label=r'FEM$_u$')

    plt.plot(T_qs_struct[:, idx], y_qs_struct,
             color="black", marker="<", markevery=5, markersize=5, linestyle="-", linewidth=2,
             label=r'FEM')

    # ------ format y axis ------ #
    # yll = y.min()
    # yhl = y.max()
    # axes.set_ylim([yll, yhl])
    # axes.set_yticks(np.linspace(yll, yhl, 5))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
    # axes.yaxis.set_major_formatter(xfmt)

    # plt.yscale('log')
    # ------ format x axis ------ #
    # plt.xlim(-0.05, 1.05)

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    plt.title(title)
    plt.xlabel(r'Temperature')
    plt.ylabel(r'$\overline{y}$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure

for cross_section, x_cut in cross_sections.items():
    # title = r'$U_x$' + f' LBM vs FEM cross-section {cross_section}'
    title = ''  # skip title for latex report
    fig_name = f'plots/ux_LBM_vs_ToolBox_{mesh}_Re{Re}_Pr{Pr}_cross_section_{cross_section}.pdf'

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    axes = plt.gca()

    plt.plot(ux_lb1[:, x_cut*1]*1, y_lb1,
             color="black", marker="", markevery=3, markersize=5, linestyle=":", linewidth=2,
             label='LBM - lattice I')

    plt.plot(ux_lb2[:, x_cut*2]*2, y_lb2,
             color="black", marker="", markevery=4, markersize=5, linestyle="-.", linewidth=2,
             label='LBM - lattice II')

    plt.plot(ux_lb3[:, x_cut*4]*4, y_lb3,
             color="black", marker="", markevery=5, markersize=5, linestyle="--", linewidth=2,
             label='LBM - lattice III')

    idx = list(cross_sections.keys()).index(cross_section)
    # plt.plot(ux_qs[:, idx], y_qs,
    #          color="black", marker=">", markevery=5, markersize=5, linestyle="-", linewidth=2,
    #          label=r'FEM$_u$')

    plt.plot(ux_qs_struct[:, idx], y_qs_struct,
             color="black", marker="<", markevery=5, markersize=5, linestyle="-", linewidth=2,
             label=r'FEM')

    # ------ format y axis ------ #
    # yll = y.min()
    # yhl = y.max()
    # axes.set_ylim([yll, yhl])
    # axes.set_yticks(np.linspace(yll, yhl, 5))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
    # axes.yaxis.set_major_formatter(xfmt)

    # plt.yscale('log')
    # ------ format x axis ------ #
    # plt.xlim(-0.05, 1.05)

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    plt.title(title)
    plt.xlabel(r'$U_{x}$')
    plt.ylabel(r'$\overline{y}$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure

print('bye')
