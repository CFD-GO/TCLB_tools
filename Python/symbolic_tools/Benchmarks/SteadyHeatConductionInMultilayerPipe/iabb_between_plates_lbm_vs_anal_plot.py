

from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, get_r_from_xy
from DataIO.helpers import calc_L2, calc_mse
import pwd
from Benchmarks.PoiseuilleFlow.pipe_poiseuille_plots import cntr_plot, slice_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot
from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import HeatConductionBetweenTwoPlates
from DataIO.VTIFile import VTIFile
import os
import pwd
import numpy as np

Heff=60.5
conductivity = 0.1
kin_visc = 0.1

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IABB_ruraWrurze')
case_folder = f'abb_ruraWrurze_Dirichlet_Cumulants_k_{conductivity}_nu_{kin_visc}_effdiam_{eff_pipe_diam}'


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

hcbp = HeatConductionBetweenTwoPlates(T0=11, T2=10, Heff=Heff)

x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)

T_anal = np.zeros((ny, nx))
r_anal = np.zeros((ny, nx))

cuttoff_y2 = Heff - 1
cuttoff_y0 = 0 + 1
for i in range(ny):
    for j in range(nx):
        T_anal[i, j] = hcbp.get_T_profile(i)
        if i < cuttoff_y0 or i > cuttoff_y2:
            T_anal[i, j] = np.nan


not_nan_mask = ~np.isnan(T_anal)
T_anal_masked = T_anal[not_nan_mask]
T_num_slice_masked = T_num_slice[not_nan_mask]
r_anal_masked = r_anal[not_nan_mask]

T_mse = calc_mse(T_anal_masked, T_num_slice_masked)
T_L2 = calc_L2(T_anal_masked, T_num_slice_masked)

# print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
# print(f"T_L2={T_L2[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")

cntr_plot(T_anal, T_num_slice, xx, yy, conductivity, Heff)
# # 2D clip
# r_anal = r_anal[:, 63]
# u_anal = u_anal[:, 63]
# uz_num_slice = uz_num_slice[:, int(nx / 2)]
# slice_plot(u_anal, uz_num_slice, kin_visc, effdiam, r_anal)

print('bye')
