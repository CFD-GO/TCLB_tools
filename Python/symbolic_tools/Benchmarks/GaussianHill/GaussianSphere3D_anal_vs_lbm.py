from DataIO.helpers import find_oldest_iteration, calc_mse, calc_L2, calc_L2_per_element
import matplotlib.pyplot as plt
import os
import numpy as np
import pwd
from DataIO.VTIFile import VTIFile
from Benchmarks.GaussianHill.GaussianHillAnal import GaussianHillAnal, prepare_anal_data_ADE_Gaussion_Hill_3D

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import plotly.graph_objects as go
import numpy as np
from vedo import *

from sympy.matrices import Matrix
import time

start = time.process_time()
# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir

ux = 0
time_SI = 1  # --> k=1e-5
# time_SI = 100  # --> k=1e-3
conductivity_SI = 4.0
Sigma02 = 100
domain_size_SI = 256
lattice_size = 256

iterations = 400000
dx = domain_size_SI / lattice_size

C0 = 1.
reference_level = 10.
X0 = Matrix([lattice_size/2., lattice_size/2., lattice_size/2.])
conductivity = time_SI*conductivity_SI/(iterations*dx*dx)
# iterations = conductivity_SI*time_SI*10000

main_dir = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianSphere_3D')
plot_dir = 'ADE_GaussianSphere3D'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

def plot_err_field(X, Y, Z, values, fig_name):
    values.min()
    values.max()
    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=values.flatten(),
        isomin=1.2*values.min(),
        # isomin=0.90 * values.max(),
        isomax=0.95 * values.max(),
        # isomin=-0.1,
        # isomax=0.8,
        opacity=0.1,  # needs to be small to see through all surfaces
        surface_count= 21,  # 21 needs to be a large number for good volume rendering
    ))
    fig.show()

def plot_err_field_vedo(X, Y, Z, scalar_field, fig_name):
    isomin = 1. * scalar_field.min()
    isomax = 1. * scalar_field.max()
    # isomin = 10
    # isomax = 11

    vol = Volume(scalar_field)
    # vol.addScalarBar3D()
    vol.addScalarBar3D(sy=1.5, title='height is the scalar')

    # generate an isosurface the volume for each thresholds
    ts = np.linspace(isomin, isomax, 6)
    # ts = [10.1, 10.25, 10.4]
    # Use c=None to use the default vtk color map. isos is of type Mesh
    isos = Volume(scalar_field).isosurface(threshold=ts).opacity(0.4)
    text1 = Text2D(f'Make a Volume from numpy.mgrid', c='blue')
    text2 = Text2D(f'its isosurface representation\nvmin={isomin:.2e}, vmax={isomax:.2e}', c='dr')

    print('numpy array from Volume:',
          vol.getPointArray().shape,
          vol.getDataArray().shape)

    # show([(vol, text1), (isos, text2)], N=2, azimuth=10, zoom=1.3)
    destination_path = os.path.join(plot_dir, fig_name)
    write(vol, destination_path)
    print(f'zapisano w: {destination_path}')


def get_t_err(main_folder, collision_type):
    case_name = f'{collision_type}_ux_{ux:.2e}_k_{conductivity:.2e}_iterations_{iterations}_sigma2_{Sigma02}_size_{lattice_size}lu'
    folder = os.path.join(main_folder, case_name)

    oldest = find_oldest_iteration(folder)
    # filename_vtk = f'sample_VTK_P00_{oldest}.vti'
    filename_vtk = f'{case_name}_VTK_P00_{oldest}.vti'

    filepath_vtk = os.path.join(folder, filename_vtk)
    vti_reader = VTIFile(filepath_vtk)
    T_num = vti_reader.get("T")

    ySIZE, xSIZE, zSIZE = T_num.shape
    assert ySIZE == xSIZE == zSIZE == lattice_size

    # dump_file_path = os.path.join(main_folder, f'dumps', f'sample.npy')

    assert lattice_size == domain_size_SI
    dump_file_path = os.path.join(main_folder, f'dumps', f'ux_{ux:.2e}_k_{conductivity:.2e}_iterations_{iterations}_sigma2_{Sigma02}_size_{lattice_size}lu.npy')
    gha = GaussianHillAnal(C0, X0, Sigma02, conductivity, Matrix([0, 0, 0]), D=3.)
    xx, yy, zz, T_anal = prepare_anal_data_ADE_Gaussion_Hill_3D(
        gha, ux, iterations,
        lattice_size, lattice_size, lattice_size, dump_file_path,
        shall_recalculate_results=False, reference_level=reference_level)

    # T_num = T_num[:, :, 127]
    # T_anal = T_anal[:, :, 0]
    T_err_field = T_anal - T_num
    T_L2 = calc_L2_per_element(T_anal, T_num)
    T_L2_sum = calc_L2(T_anal, T_num)
    T_mse_sum = calc_mse(T_anal, T_num)
    # print(f"T_mse={T_mse[g]:.2e} for grid {xSIZE} x {xSIZE} [lu]")
    # print(f"{collision_type} T_L2={T_L2[g]:.5e} for k = {conductivities[g]}")

    print("------------------------------------ PLOT err field------------------------------------")
    # half_size = int(lattice_size/2)
    # T_slice = T_err_field[:half_size, :half_size, :half_size]
    plot_err_field_vedo(xx, yy, zz, T_err_field, fig_name=f'{case_name}_err_field.vti')
    plot_err_field_vedo(xx, yy, zz, T_anal, fig_name=f'{case_name}_anal_field.vti')
    plot_err_field_vedo(xx, yy, zz, T_num, fig_name=f'{case_name}_num_field.vti')
    # plot_err_field(xx, yy, zz, T_err_field, fig_name)
    return T_L2
    # return T_mse


T_err_L2_BGK = get_t_err(main_dir, 'BGK')
T_err_L2_CM = get_t_err(main_dir, 'CM')
T_err_L2_CM_HIGHER = get_t_err(main_dir, 'CM_HIGHER')
T_err_L2_Cumulants = get_t_err(main_dir, 'Cumulants')

print(f'\n\n Done in {time.process_time() - start} [s].')
