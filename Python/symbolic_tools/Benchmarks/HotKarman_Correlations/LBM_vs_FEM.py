
from DataIO.VTIFile import VTIFile
import os
import pandas as pd
import pwd
import numpy as np
import matplotlib.pyplot as plt

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotBarman3D')
sizer = 2
# readme:
# sizer 1 = LB mesh: 1000x150x3 -- sample 150 y points from toolbox
# sizer 2 = LB mesh: 2000x300x6 -- sample 300 y points from toolbox
# sizer 4 = LB mesh: 4000x600x12 -- sample 600 y points from toolbox

filename_vtk = f'HotKarman3D_template_sizer_{sizer}_nu_0.03_k_0.003_VTK_P00_00600000.pvti'
filepath_pvti = os.path.join(main_folder, f'lbm_template_sizer_{sizer}_nu_0.03_k_0.003', filename_vtk)

vti_reader = VTIFile(filepath_pvti, parallel=True)

T_num = vti_reader.get("T")
ny, nx, nz = T_num.shape

T_num_slice = T_num[:, :, 1]
# # convert to pandas format
# data = T_num_slice
# # pdT = pd.DataFrame(data=data[1:, 1:],  # values
# #                    index=data[1:, 0],  # 1st column as index
# #                    columns=data[0, 1:])  # 1st row as the column names
#
# pdT2 = pd.DataFrame(T_num_slice)


# Read data from toolbox
filepath_csv = os.path.join(main_folder, 'from_toolbox', f'Temperature_{int(150*sizer)}_results_clmax_3_Re_10_Pr_10.csv')
# pdT = pd.read_csv(filepath_csv, delimiter=",", header=None)
pdT2 = pd.read_csv(filepath_csv, delimiter=",")
# T_toolbox = np.genfromtxt(filepath_csv, delimiter=',', names=True)
T_toolbox = np.loadtxt(filepath_csv, delimiter=',', skiprows=1)
print("data reading complete, lets plot it!")

y = np.linspace(start=0, stop=ny, num=ny, endpoint=False)


def make_plot(y, T_lbm, T_qs_toolbox, title, fig_name):
    # x = T_toolbox[:, 0]

    # -------------------- make dummy plot --------------------

    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(14, 8))

    axes = plt.gca()
    plt.plot(T_lbm, y,
             color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
             label='LBM')

    plt.plot(T_qs_toolbox, y,
             color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=2,
             label='toolbox')

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
    plt.xlabel(r'$Temperature$')
    plt.ylabel(r'$y_{coord}$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure


# title = f'LBM_vs_ToolBox.png'
# make_plot(y, T_num_slice[:, 270], T_toolbox[:, 0], title)

x_cuts = np.array([270, 300, 330, 960, 999])*sizer

for i in range(len(x_cuts)):
    title = f'LBM vs ToolBox \n {int(1000*sizer)}x{int(150*sizer)}[lu] x_cut={x_cuts[i]}'
    fig_name = f'LBM_vs_ToolBox_size={int(1000*sizer)}x{int(150*sizer)}[lu]_x_cut={x_cuts[i]}'
    make_plot(y, T_num_slice[:, x_cuts[i]], T_toolbox[:, i], title, fig_name)
