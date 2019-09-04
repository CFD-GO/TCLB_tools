import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pwd
from DataIO.VTIFile import VTIFile
# from Benchmarks.HotKarman_Correlations.LBM_vs_FEM_v2 import read_data_from_LBM

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotBarman3D')
clmax = 1

Re = 10
Pr = 1000
def read_data_from_LBM_for_cntr_plot(size):
    # filename_vtk = f'HotKarman3D_template_sizer_{sizer}_nu_0.03_k_0.003_VTK_P00_01000000.pvti'
    filename_vtk = f'HotKarman3D_template_sizer_{size}_Re{Re}_Pr{Pr}_VTK_P00_00980000.pvti'
    filepath_pvti = os.path.join(main_folder, f'keepU_sizer_{size}_Re{Re}_Pr{Pr}', filename_vtk)

    vti_reader = VTIFile(filepath_pvti, parallel=True)

    T_num = vti_reader.get("T")
    [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    ny, nx, nz = T_num.shape

    ux_num_slice = ux_num[:, :, 1]
    T_num_slice = T_num[:, :, 1]

    x_grid = np.linspace(start=0, stop=nx, num=nx, endpoint=False)
    y_grid = np.linspace(start=0, stop=ny, num=ny, endpoint=False)
    z_grid = np.linspace(start=0, stop=nz, num=nz, endpoint=False)

    xx, yy = np.meshgrid(x_grid, y_grid)

    return T_num_slice, ux_num_slice, xx, yy

T_lb_slice1, ux_lb1, xx_lb1, yy_lb1 = read_data_from_LBM_for_cntr_plot(size=4)
# T_lb_slice2, ux_lb2, y_lb2 = read_data_from_LBM(size=2)
# T_lb_slice3, ux_lb3, y_lb3 = read_data_from_LBM(size=4)

# ny, nx, nz = T_lb_slice1.shape
# ny, nx, nz = T_lb_slice2.shape
# ny, nx, nz = T_lb_slice3.shape


print("---------- PLOTTING -------------")

# T_num_slice = T_num[:, :, 1]

fig = plt.figure(figsize=(12, 8))
# ax = fig.gca(projection='3d')
ax = fig.gca()
# alpha=1, rstride=1, cstride=1)
# ax.plot_surface(xx, yy, T_err_field, cmap='winter', linewidth=0.5, antialiased=True, zorder=0.5, label='T_err_field')
# ax.plot_surface(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, zorder=0.25, label='T_num')
# ax.plot_surface(xx, yy, T_anal,  cmap='autumn', linewidth=0.5, antialiased=True, zorder=0.1, label='T_anal')
# ax.contourf(xx, yy, T_num_slice,  cmap='summer', linewidth=0.5, antialiased=True, label='T_num')

# cntr = ax.pcolormesh(xx_lb1, yy_lb1, T_lb_slice1, cmap='coolwarm', label='T_num')  # this one has smooth colors
# cntr = ax.contourf(xx_lb1, yy_lb1, T_lb_slice1, cmap='coolwarm', antialiased=True)  # this one is has step colors

cntr = ax.pcolormesh(xx_lb1, yy_lb1, ux_lb1, cmap='coolwarm', label='T_num')  # this one has smooth colors

ax.set_xlabel('X')
ax.set_ylabel('Y')
# ax.set_zlabel('Z')
ax.set_aspect('equal')

# Customize the z axis.
# ax.set_zlim(-.1, 1.05)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.1e'))


# Add a color bar which maps values to colors.

fig.colorbar(cntr, shrink=0.5, aspect=5)

plt.title(f'Laplace benchmark\n ')
# plt.grid(True)

# fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
# ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)
fig_name = f'Plot_from_3D_data.png'

fig.savefig(fig_name, bbox_inches='tight')
plt.show()


print('bye')
