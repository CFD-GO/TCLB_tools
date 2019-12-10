

from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, get_r_from_xy
from DataIO.helpers import calc_L2, calc_mse
import pwd
from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot, slice_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import HeatConductionBetweenTwoPlates
from DataIO.VTIFile import VTIFile
import os
import pwd
import numpy as np

H = 63# [15, 31, 63]
q = 0.5  # [0.25, 0.5, 0.75]
expected_wall_location = 1.5 - q  # location of bottom wall, where h = 0
Heff = H - 2 * expected_wall_location

conductivity = 0.1
collision_type = 'CM_HIGHER'
bb_type = 'abb'

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IABB_plates')
case_folder = strip_folder_name(f'{bb_type}_plates_Dirichlet_{collision_type}_k_{conductivity}_effdiam_{Heff}')


def read_data_from_lbm(_case_folder):
    oldest = find_oldest_iteration(_case_folder)
    filename_vtk = get_vti_from_iteration(_case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(_case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    T_num = vti_reader.get("T")
    T_num_slice = T_num[:, :, 1]

    # [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    # ny, nx, nz = T_num.shape
    # uz_num_slice = uz_num[:, :, 1]
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return T_num_slice


T_num_slice = read_data_from_lbm(os.path.join(main_folder, case_folder))
ny, nx = T_num_slice.shape

# -------- anal solution ---------------

hcbp = HeatConductionBetweenTwoPlates(T0=11, T2=10, Heff=Heff, y0=expected_wall_location)

x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
# np.insert(y_grid, len(y_grid)-2, Heff, axis=0)

xx, yy = np.meshgrid(x_grid, y_grid)

T_anal = np.zeros((ny, nx))
# r_anal = np.zeros((ny, nx))

cuttoff_y2 = Heff - 1
cuttoff_y0 = 0 + 1
for i in range(ny):
    for j in range(nx):
        T_anal[i, j] = hcbp.get_T_profile(yy[i][j])
        if i < cuttoff_y0 or i > cuttoff_y2:
            T_anal[i, j] = np.nan


not_nan_mask = ~np.isnan(T_anal)
T_anal_masked = T_anal[not_nan_mask]
T_num_slice_masked = T_num_slice[not_nan_mask]
yy_masked = yy[not_nan_mask]

T_mse = calc_mse(T_anal_masked, T_num_slice_masked)
T_L2 = calc_L2(T_anal_masked, T_num_slice_masked)

# cntr_plot(T_anal, T_num_slice, xx, yy, conductivity, Heff, title=case_folder)

# # 2D clipy_grid
# yy_masked = yy_masked[:,  int(nx / 2)]
T_anal_slice = T_anal[:,  int(nx / 2)]
T_num_slice = T_num_slice[:, int(nx / 2)]
not_nan_mask = ~np.isnan(T_anal_slice)
T_anal_masked = T_anal_slice[not_nan_mask]
T_num_slice_masked = T_num_slice[not_nan_mask]
y_masked = y_grid[not_nan_mask]

slice_plot(T_anal_masked, T_num_slice_masked, y_masked, title=f'{case_folder}_q{q}')


print(f"T_mse={T_mse:.2e} for k{conductivity}_q{q}_effdiam_{Heff}")
print(f"T_L2={T_L2:.2e} for k{conductivity}_q{q}_effdiam_{Heff}")

print('bye')
