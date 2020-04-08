import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.VTIFile import VTIFile
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration
import os
import pandas as pd
import numpy as np
import pwd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import re
import warnings
from fractions import Fraction

home = pwd.getpwuid(os.getuid()).pw_dir
local_logs_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_square_advection')

#
# wd = os.getcwd()
# wd = os.path.dirname(wd)  # go level up

def make_plot(data, plot_name):

    # filename = 'HotKarman3D_template_sizer_1_nu_0.03_k_0.003_VTK_P00_00600000.pvti'
    # filepath = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', '3D_multiGPU', filename)


    # vti_reader = VTIFile(file_path, parallel=False)
    # T_num = vti_reader.get("T")

    T_num = data
    ny, nx, nz = T_num.shape

    # U = vti_reader.get("U", vector=True)

    # convert to pandas format
    # data = T_num
    # pdT = pd.DataFrame(data=data[1:, 1:],  # values
    #                    index=data[1:, 0],  # 1st column as index
    #                    columns=data[0, 1:])  # 1st row as the column names
    #
    # pdT2 = pd.DataFrame(data)


    print("---------- PLOTTING -------------")
    x_grid = np.linspace(start=0, stop=nx, num=nx, endpoint=False)
    y_grid = np.linspace(start=0, stop=ny, num=ny, endpoint=False)
    z_grid = np.linspace(start=0, stop=nz, num=nz, endpoint=False)
    xx, yy = np.meshgrid(x_grid, y_grid)
    T_num_slice = T_num[:, :, 1]

    fig = plt.figure(figsize=(12, 8))
    # ax = fig.gca(projection='3d')
    ax = fig.gca()
    # alpha=1, rstride=1, cstride=1)
    # ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5, label='T_err_field')
    # ax.plot_surface(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, label='T_num')
    # ax.plot_surface(xx, yy, T_anal,  cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.1, label='T_anal')
    # ax.contourf(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, label='T_num')

    # gist_rainbow, brg, hsv, seismic, RdBu
    cntr = ax.pcolormesh(xx, yy, T_num_slice, cmap='Greys', label='T_num', vmin=9.975,
                         vmax=11.025)  # this one has smooth colors

    # cntr = ax.pcolormesh(xx, yy, T_num_slice, cmap='twilight', label='T_num', vmin=9.95, vmax=11.05)  # this one has smooth colors
    # cntr = ax.contourf(xx, yy, T_num_slice, cmap='coolwarm', antialiased=True)  # this one is has step colors

    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    # Customize the z axis.
    # ax.set_zlim(-.1, 1.05)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter(FormatStrFormatter('%.1e'))


    # Add a color bar which maps values to colors.
    # fig.colorbar(cntr, shrink=0.5, aspect=5)

    # plt.title(f'Laplace benchmark\n ')

    # # Major ticks every 20, minor ticks every 5
    # major_ticks = np.arange(0, nx, nx/5)
    # minor_ticks = np.arange(0, nx, nx/10)
    #
    # ax.set_xticks(major_ticks)
    # ax.set_xticks(minor_ticks, minor=True)
    # # ax.set_yticks(major_ticks)
    # # ax.set_yticks(minor_ticks, minor=True)
    #
    # # And a corresponding grid
    # ax.grid(which='both')
    #
    # # Or if you want different settings for the grids:
    # ax.grid(which='minor', alpha=0.2)
    # ax.grid(which='major', alpha=0.5)
    # plt.grid(True)  # or use default grid

    # fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
    # ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)

    match = re.search(r'k_(0.\d{1,4})|k_(1e-\d{2})', plot_name, re.IGNORECASE)
    k_dir = match.group(0)

    subdir = os.path.join('plots', f'{k_dir}')
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    fig_name = f'{subdir}/{plot_name}.png'

    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()
    plt.close(fig)  # close the figure



def read_data_from_LBM(case_folder):
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    # match = re.search(r'Pr_?(\d{1,4})_', file, re.IGNORECASE)
    # Pr = float(match.group(1))
    T_num = vti_reader.get("T")

    name_for_plot = re.sub(r"VTK_P00_", '', filename_vtk)
    name_for_plot = re.sub(r".pvti", '', filename_vtk)
    return T_num, name_for_plot
    # [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    # ny, nx, nz = T_num.shape
    #
    # ux_num_slice = ux_num[:, :, 1]
    # T_num_slice = T_num[:, :, 1]
    #
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    # return T_num_slice, ux_num_slice, y


for root, dirs, files in os.walk(local_logs_folder):
    for directory in dirs:
        case_folder = os.path.join(root, directory)
        data, name_for_plot = read_data_from_LBM(case_folder)
        make_plot(data, plot_name=name_for_plot)

    # match = re.search(r'Pr_?(\d{1,4})_', file, re.IGNORECASE)
    # Pr = float(match.group(1))
    # for file in files:
    #     if file.endswith('.pvti') and 'toolbox' not in root:
    #         print(file)


print('bye')
