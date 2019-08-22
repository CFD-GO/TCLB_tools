import matplotlib.pyplot as plt
import numpy as np
import os
from Benchmarks.PoiseuilleFlow.PoiseuilleAnal import calc_gx_in_pipe, OnePhasePoiseuilleAnalInPipe
from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, strip_folder_name, eat_dots_for_texmaker, get_r_from_xy
from DataIO.helpers import calc_L2, calc_mse
import pwd
from DataIO.VTIFile import VTIFile
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import os
from DataIO.VTIFile import VTIFile
from Benchmarks.PoiseuilleFlow.pipe_poiseuille_plots import cntr_plot, slice_plot

rho = 1.
kin_visc = 0.1
mu = rho * kin_visc

effdiam = 121
uc = 0.01

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_IBB_PoiseuillePipe')

case_folder = f'ibb_PoiseuillePipe_CM_HIGHER_nu_{kin_visc}_effdiam_{effdiam}'


def read_data_from_LBM(case_folder):
    oldest = find_oldest_iteration(case_folder)
    filename_vtk = get_vti_from_iteration(case_folder, oldest, extension='.pvti')
    filepath_vtk = os.path.join(case_folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk, parallel=True)
    T_num = vti_reader.get("T")
    [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    ny, nx, nz = T_num.shape

    uz_num_slice = uz_num[:, :, 1]
    # T_num_slice = T_num[:, :, 1]
    #
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return uz_num_slice


uz_num_slice = read_data_from_LBM(os.path.join(main_folder, case_folder))
ny, nx = uz_num_slice.shape

# -------- anal solution ---------------
# y_anal = np.arange(q, effdiam, 1)
# y_anal = np.concatenate(([0], y_anal, [effdiam]))
# r_anal = np.arange(0, effdiam/2, 1)
gx = calc_gx_in_pipe(uc, rho, kin_visc, D=effdiam)
poiseuilleAnal = OnePhasePoiseuilleAnalInPipe(gx=gx, nu=kin_visc, D=effdiam)
# u_anal = np.array([poiseuilleAnal.get_u_profile(r_anal[i]) for i in range(len(r_anal))])


x0 = 63.5
y0 = 63.5

x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)

u_anal = np.zeros((ny, nx))
r_anal = np.zeros((ny, nx))
cuttoff = effdiam / 2. - 1


for i in range(ny):
    for j in range(nx):
        r = get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
        r_anal[i, j] = r
        u_anal[i, j] = poiseuilleAnal.get_u_profile(r)
        if r > cuttoff:
            u_anal[i, j] = np.nan


not_nan_mask = ~np.isnan(u_anal)
u_anal_masked = u_anal[not_nan_mask]
uz_num_slice_masked = uz_num_slice[not_nan_mask]
r_anal_masked = r_anal[not_nan_mask]

uz_mse = calc_mse(u_anal_masked, uz_num_slice_masked)
uz_L2 = calc_L2(u_anal_masked, uz_num_slice_masked)

# print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
# print(f"T_L2={T_L2[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")


cntr_plot(u_anal, uz_num_slice, kin_visc, effdiam, nx, xx, yy)
# # 2D clip
# r_anal = r_anal[:, 63]
# u_anal = u_anal[:, 63]
# uz_num_slice = uz_num_slice[:, int(nx / 2)]
# slice_plot(u_anal, uz_num_slice, kin_visc, effdiam, r_anal)

print('bye')
