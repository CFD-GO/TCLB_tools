import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.VTIFile import VTIFile


wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

filename = 'HotKarman3D_template_sizer_1_nu_0.03_k_0.003_VTK_P00_00600000.pvti'
filepath = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', '3D_multiGPU', filename)


vti_reader = VTIFile(filepath, parallel=True)
T_num = vti_reader.get("T")

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

# cntr = ax.pcolormesh(xx, yy, T_num_slice, cmap='coolwarm', label='T_num')  # this one has smooth colors
cntr = ax.contourf(xx, yy, T_num_slice, cmap='coolwarm', antialiased=True)  # this one is has step colors

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
plt.grid(True)  # or use default grid



# fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
# ax.legend([fake2Dline], [r'$Err_{T} = T_{anal} - T_{num}$'], numpoints=1)

if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/Plot_from_3D_data.png'

fig.savefig(fig_name, bbox_inches='tight')
plt.show()
# plt.close(fig)  # close the figure

print('bye')
