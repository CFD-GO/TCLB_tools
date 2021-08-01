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
from vedo import datadir, show, Volume

from sympy.matrices import Matrix
import time

start = time.process_time()
# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up
home = pwd.getpwuid(os.getuid()).pw_dir

# ux = 0
# time_SI = 1  # --> k=1e-5
# # time_SI = 100  # --> k=1e-3
# conductivity_SI = 4.0
# Sigma02 = 100
# domain_size_SI = 256
# lattice_size = 256
#
# iterations = 400000
# dx = domain_size_SI / lattice_size
#
# C0 = 1.
# reference_level = 10.
# X0 = Matrix([lattice_size/2., lattice_size/2., lattice_size/2.])
# conductivity = time_SI*conductivity_SI/(iterations*dx*dx)
# # iterations = conductivity_SI*time_SI*10000

main_dir = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_GaussianSphere_3D', 'err_fields')
plot_dir = 'ADE_GaussianSphere3D_plots'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

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


def plot_isos_vedo(vti_path, output_path):
    # generate an isosurface the volume for each thresholds

    vol = load(vti_path)
    pa = vol.getPointArray()
    da = vol.getDataArray()
    # isos is of type Mesh
    # isos = Volume(datadir+'quadric.vti').isosurface(thresholds)
    isomin = 1. * pa.min()
    isomax = 1. * pa.max()
    ts = np.linspace(isomin, isomax, 6)
    thresholds = [1E-3, 1E-4, 1E-5]
    isos = vol.isosurface(thresholds)
    show(isos,
         "click mesh and press shift-X",
         axes=1, interactive=False)

    # write(isos, output_path)
    # screenshot('output_path')
    print(f'saved in {output_path}')


for root, dirs, files in os.walk(main_dir):
    for file in files:
        # if file.endswith('.vti'):
        if file.endswith('CM_ux_0.00e+00_k_1.00e-03_iterations_400000_sigma2_100_size_256lu_err_field.vti'):
            filepath = os.path.join(root, file)
            plot_isos_vedo(filepath,
                           os.path.join(plot_dir, f'contours_bee.png')
                           )

print(f'\n\n Done in {time.process_time() - start} [s].')
