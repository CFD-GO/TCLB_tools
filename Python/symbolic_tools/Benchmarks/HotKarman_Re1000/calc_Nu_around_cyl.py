

from DataIO.helpers import find_oldest_iteration, get_vti_from_iteration, get_r_from_xy
from DataIO.helpers import calc_L2, calc_mse
# from Benchmarks.SteadyHeatConductionInMultilayerPipe.contour_and_slice_plot import cntr_plot
# from Benchmarks.SteadyHeatConductionInMultilayerPipe.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeDirichlet
from DataIO.VTIFile import VTIFile
import re
import os
import pwd
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp2d


from Benchmarks.HotKarman_Re1000.contour_and_slice_plot import cntr_plot, Nu_phi_plot

# eff_pipe_diam = 66
# eff_cyl_diam = 33
# conductivity = 0.1
# kin_visc = 0.1

Pr = 0.71
dr = 1.5
whichT = 'T'  # 'T', 'averageT'  # TODO: add mask T=11 for averageT
solver = 'TCLB'  # 'TCLB' or 'walberla

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir

if solver == 'walberla':
    main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'walberla_IABB_ruraWrurze')
    # case_folder = os.path.join('vtk_test', 'thermal_field')
    # case_folder = os.path.join('vtk_test_old', 'thermal_field')
    case_folder = os.path.join(f'vtk_eff_pipe_diam_{eff_pipe_diam}', 'thermal_field')
elif solver == 'TCLB':
    main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotKarman_Re1000')
    # case_folder = f'abb_ruraWrurze_Dirichlet_Cumulants_k_{conductivity}_nu_{kin_visc}_effdiam_{eff_pipe_diam}'
    # case_folder = f'test_batch_HotKarman_Re1000_D0_1.28e+02_nu_0.0055_U_4.30e-02_1st_order_bc_CM_HIGHER'
    case_folder = 'batch_HotKarman_Re1000_D0_6.40e+01_nu_0.003_U_4.69e-02_1st_order_bc_CM_HIGHER'
else:
    raise Exception("Choose solver [\'TCLB\' or \'walberla\'] ")


def read_data_from_lbm(_case_folder):
    if solver == 'walberla':
        oldest = find_oldest_iteration(_case_folder, extension='.vti')
        filename_vtk = get_vti_from_iteration(_case_folder, oldest, extension='.vti', prefix='simulation_step_')  # walberla
        filepath_vtk = os.path.join(_case_folder, filename_vtk)
        vti_reader = VTIFile(filepath_vtk, parallel=False)  # walberla
    elif solver == 'TCLB':
        oldest = find_oldest_iteration(_case_folder)  # TCLB
        filename_vtk = get_vti_from_iteration(_case_folder, oldest, extension='.pvti') # TCLB
        filepath_vtk = os.path.join(_case_folder, filename_vtk)
        vti_reader = VTIFile(filepath_vtk, parallel=True)  # TCLB


    else:
        raise Exception("Choose solver [\'TCLB\' or \'walberla\'] ")

    T_num = vti_reader.get(whichT)
    T_num_slice = T_num[:, :, 1]

    # [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)
    # ny, nx, nz = T_num.shape
    # uz_num_slice = uz_num[:, :, 1]
    # y = np.linspace(start=0, stop=1, num=ny, endpoint=False)
    return T_num_slice


def get_sim_params_from_folder_name(_case_folder):
    match = re.search(r'Re(\d{1,4})_', _case_folder, re.IGNORECASE)
    Re = float(match.group(1))

    match = re.search(r'_U_(\d\.\d\de\-\d\d)_', _case_folder, re.IGNORECASE)
    u = float(match.group(1))

    match = re.search(r'_D0_(\d\.\d\de\+\d\d)_', _case_folder, re.IGNORECASE)
    d = float(match.group(1))

    match = re.search(r'_nu_(0\.\d{1,4})_', _case_folder, re.IGNORECASE)
    kin_visc = float(match.group(1))

    np.testing.assert_allclose(Re, u*d/kin_visc, rtol=1e-1)
    return u, d, kin_visc


T_num_slice = read_data_from_lbm(os.path.join(main_folder, case_folder))
ny, nx = T_num_slice.shape

u, d, kin_visc = get_sim_params_from_folder_name(os.path.join(main_folder, case_folder))
conductivity = kin_visc/Pr
r = d/2.  # true for I(A)BB

x0 = 0.25 * nx  # x-center of the cylinder/pipe
y0 = 0.5 * ny   # y-center of the cylinder/pipe

# generate points where data needs to be interpolated in cartesian coord
x_grid = np.linspace(0, nx, nx, endpoint=False) + 0.5
y_grid = np.linspace(0, ny, ny, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)

# transform them to radial coord
phi_deg = np.linspace(0, 360, num=120)
x_r = np.array(x0 + (r+dr)*np.cos(np.deg2rad(phi_deg)))
y_r = np.array(y0 + (r+dr)*np.sin(np.deg2rad(phi_deg)))

# https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy

# interpolation 1 - poor results
# interpol_fun = interp2d(x_grid, y_grid, T_num_slice, kind='cubic')
# Tr = interpol_fun(x_r, y_r)
# Tr = Tr.diagonal()
# Nu_phi_plot(Tr, phi_deg, title=r'$T(\phi)$', file_name='interp2d')

# interpolation 2
ipoints_xr, ipoints_yr = np.meshgrid(x_r, y_r)  # points where data needs to be interpolated
known_data = np.array([xx.flatten(), yy.flatten()]).T
Tr = griddata(known_data, T_num_slice.flatten(), (ipoints_xr, ipoints_yr), method='linear')
Tr = Tr.diagonal()
Nu_phi_plot(Tr, phi_deg, title=r'$T(\phi)$', file_name=f'T_dr{dr}_{whichT}')

# T_inf = 10
# q = conductivity*(Tr - T_inf)
# L = 1  # we are already on slice
# S = np.pi * d * L
# h = q/(S*(Tr-T_inf))
# Nu = h*d/conductivity


def calc_Nu(q_conv, k, D, L):
    T_surf = 11
    T_inf = 10
    Surface = np.pi * D * L
    ht_coeff_experimental = q_conv / (Surface * (T_surf - T_inf))
    Nu = ht_coeff_experimental * D / k
    Nu_jfm = q_conv*D/(k*L*(T_surf - T_inf))  # TODO:
    return Nu, Nu_jfm

T_inf = 10
q = conductivity*(Tr - T_inf)
Nu, Nu_jfm = calc_Nu(q, conductivity, d, 1)

# TODO: policz calke z Nu i porownaj z logiem
Nu_phi_plot(Nu_jfm, phi_deg, title=r'$Nu(\phi)$', file_name=f'Nu_dr{dr}_{whichT}')


# x_grid2 = np.linspace(0, nx, 0.2*nx, endpoint=False) + 0.5
# y_grid2 = np.linspace(0, ny, 0.2*ny, endpoint=False) + 0.5
# Z2 = interpol_fun(x_grid2, y_grid2)
# cntr_plot(None, Z2, x_grid2, y_grid2)

# T_anal = np.zeros((ny, nx))
# r_anal = np.zeros((ny, nx))
#
# for i in range(ny):
#     for j in range(nx):
#         r = get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
#         r_anal[i, j] = r
#         # T_anal[i, j] = pwp.get_temperature_r(r)
#         # T_anal[i, j] = 10 # TODO
#         if r < cuttoff_r0:
#             T_anal[i, j] = np.nan
#         # if r < cuttoff_r0 or r > cuttoff_r2:
#         #     T_anal[i, j] = np.nan
#
#
# not_nan_mask = ~np.isnan(T_anal)
# T_anal_masked = T_anal[not_nan_mask]
# T_num_slice_masked = T_num_slice[not_nan_mask]
# r_anal_masked = r_anal[not_nan_mask]


# show the intepolation points
points_be_shown_on_plot = np.concatenate(([x_r], [y_r]), axis=0)
cntr_plot(None, T_num_slice, xx, yy, title='interpolation_points_for_T', extra_points=points_be_shown_on_plot)

print('bye')
