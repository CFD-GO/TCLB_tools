# 2D, steady state case:
# Rectangular plate, temperature on one of the walls is given by a quadratic function.


# References:
# https://www.wire.tu-bs.de/lehre/ws15/pde1/lecture_2.pdf
# https://www.wire.tu-bs.de/lehre/ws15/pde1/lecture_3.pdf

from joblib import Parallel, delayed
import multiprocessing
from typing import NamedTuple
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
import sympy as sp
from sympy.abc import x, y, t
from sympy import sin, pi, integrate
from typing import Callable
from dataclasses import dataclass
import time

import os
import pwd


@dataclass
class InputForLaplace2DAnalytical:
    x_low: float
    x_high: float
    y_low: float
    y_high: float
    step: float
    my_fun: Callable
    n_fourier_terms: int


def analytical_laplace_2d(input_config: InputForLaplace2DAnalytical):
    lim = (x, input_config.x_low, input_config.x_high)
    Lx = input_config.x_high - input_config.x_low
    Ly = input_config.y_high - input_config.y_low
    u_sol = 0
    c = [0]

    start = time.process_time()
    print("---------- Calculating Fourier coeff -------------")
    for k in range(1, input_config.n_fourier_terms):  # skip zero
        result = 2 / (Lx * sp.sinh(k * pi * Ly / Lx)) * integrate(input_config.my_fun * sin(k * pi * x / Lx), lim)
        c.append(result)
        u_sol += c[k] * sp.sinh(k * pi * y / Lx) * sp.sin(k * pi * x / Lx)

    print("---------- Calculating Values -------------")

    nx = int((input_config.x_high + input_config.x_low) / input_config.step)
    ny = int((input_config.y_high + input_config.y_low) / input_config.step)
    x_grid = np.linspace(input_config.x_low, input_config.x_high, nx, endpoint=True)
    y_grid = np.linspace(input_config.y_low, input_config.y_high, ny, endpoint=True)
    # x_grid = np.linspace(input_config.x_low, input_config.x_high, nx, endpoint=False) + 0.5
    # y_grid = np.linspace(input_config.y_low, input_config.y_high, ny, endpoint=False) + 0.5
    xx, yy = np.meshgrid(x_grid, y_grid)
    zz = np.zeros((ny, nx))

    # -------------------------------------
    # def process_input(i):
    #     for j in range(nx):
    #         # Z[i][j] = my_fun.subs({'x': X[i][j]})
    #         Z[i][j] = u_sol.subs({'x': X[i][j], 'y': Y[i][j]})
    #

    # num_cores = multiprocessing.cpu_count()
    # print(num_cores)
    # Parallel(n_jobs=num_cores, require='sharedmem')(delayed(process_input)(i) for i in range(ny))
    # -------------------------------------

    # Row-major order is also known as the C order, as the C programming language uses it.
    # New NumPy arrays are by default in row-major order.
    # my_hacked_fun = -4 * x * (x - input_config.x_high)/(input_config.x_high*input_config.x_high)
    for i in range(ny):
        print(f"=== Calculating i/ny: {i}/{ny}  ===")
        for j in range(nx):
            # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
            # zz[i][j] = my_hacked_fun.subs({'x': xx[i][j]})
            zz[i][j] = u_sol.subs({'x': xx[i][j], 'y': yy[i][j]})

    print(f'\n\n Done in {time.process_time() - start} [s].')
    return xx, yy, zz


def prepare_anal_data_new(ySIZE, xSIZE, folder):
    # two functions are needed - one for ABB(with zeros in 1 and 63) and second for EQ (zeros in 0.5 and 63.5)
    def analytical_laplace_2d_new(is_bc_wet_node_type, n_fourier_terms):
        x1 = 0
        y1 = 0
        y2 = ySIZE
        ym = 1
        if is_bc_wet_node_type:
            x2 = xSIZE - 1
        else:
            x2 = xSIZE - 2

        xm = 0.5 * (x1 + x2)
        a = ym / ((xm - x1) * (xm - x2))
        my_fun = a * (x - x1) * (x - x2)

        lim = (x, x1, x2)
        Lx = x2 - x1
        Ly = y2 - y1
        u_sol = 0
        c = [0]

        start = time.process_time()
        print("---------- Calculating Fourier coeff -------------")
        for k in range(1, n_fourier_terms):  # skip zero
            result = 2 / (Lx * sp.sinh(k * pi * Ly / Lx)) * integrate(my_fun * sin(k * pi * x / Lx), lim)
            c.append(result)
            u_sol += c[k] * sp.sinh(k * pi * y / Lx) * sp.sin(k * pi * x / Lx)

        print("---------- Calculating Values -------------")

        nx = x2
        ny = ySIZE
        x_grid = np.linspace(start=x1, stop=xSIZE, num=nx, endpoint=True)
        y_grid = np.linspace(start=y1, stop=xSIZE, num=ny, endpoint=True)
        xx, yy = np.meshgrid(x_grid, y_grid)
        zz = np.zeros((ny, nx))

        # Row-major order is also known as the C order, as the C programming language uses it.
        # New NumPy arrays are by default in row-major order.
        # my_hacked_fun = -4 * x * (x - input_config.x_high)/(input_config.x_high*input_config.x_high)
        for i in range(ny):
            print(f"=== Calculating i/ny: {i}/{ny}  ===")
            for j in range(nx):
                # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
                zz[i][j] = u_sol.subs({'x': xx[i][j], 'y': yy[i][j]})

        print(f'\n\n Done in {time.process_time() - start} [s].')
        return xx, yy, zz

    n_fourier = 2
    is_bc_wet_node_type = None

    upper_folder = os.path.dirname(folder)
    dump_fname = None
    if 'abb' in folder:
        is_bc_wet_node_type = False
        dump_fname = os.path.join(upper_folder, f'n_fourier{n_fourier}', f'T_anal_abb_x{xSIZE}y{ySIZE}.npy')
    else:
        is_bc_wet_node_type = True
        dump_fname = os.path.join(upper_folder, f'n_fourier{n_fourier}', f'T_anal_eq_x{xSIZE}y{ySIZE}.npy')

    if os.path.isfile(dump_fname):
        print(f'{dump_fname} found, loading results from disc')
        (xx, yy, T_anal) = np.load(dump_fname)
        # T_anal = np.load(dump_fname)
    else:
        print(f'{dump_fname} not found, starting calculations')
        xx, yy, T_anal = analytical_laplace_2d_new(is_bc_wet_node_type, n_fourier)
        fourier_folder = os.path.dirname(dump_fname)
        if not os.path.exists(fourier_folder):
            os.makedirs(fourier_folder)
        np.save(dump_fname, (xx, yy, T_anal))

    return xx, yy, T_anal


def prepare_anal_data(ySIZE, xSIZE, folder):
    x1 = 0.5
    x2 = xSIZE - 0.5
    xm = 0.5 * (x1 + x2)
    ym = 1
    if 'abb' in folder:
        a = ym / ((xm - (x1+0.5)) * (xm - (x2-0.5)))
        my_fun = a * (x - (x1+0.5)) * (x - (x2-0.5))
    else:
        a = ym / ((xm - x1) * (xm - x2))
        my_fun = a * (x - x1) * (x - x2)

    n_fourier = 2
    anal_input = InputForLaplace2DAnalytical(
        x_low=x1, x_high=x2,
        y_low=x1, y_high=x2,
        step=1, my_fun=my_fun, n_fourier_terms=n_fourier
    )

    is_bc_wet_node_type = None
    if 'abb' in folder:
        is_bc_wet_node_type = False
    else:
        is_bc_wet_node_type = True


    # # SIN fun version
    # A = sp.pi / (x2 - x1)
    # B = x1
    # my_sin_fun = sp.sin(A * (x - B))
    #
    # n_fourier = 4
    # anal_input = InputForLaplace2DAnalytical(
    #     x_low=x1, x_high=x2,
    #     y_low=x1, y_high=x2,
    #     step=1, my_fun=my_sin_fun, n_fourier_terms=n_fourier
    # )

    # dump_fname = os.path.join(main_folder, f'n_fourier{n_fourier}', f'T_anal_x{xSIZE}y{ySIZE}.npy')  # old

    # two functions are needed - one for ABB(with zeros in 1 and 63) and second for EQ (zeros in 0.5 and 63.5)
    upper_folder = os.path.dirname(folder)
    dump_fname = None
    if 'abb' in folder:
        dump_fname = os.path.join(upper_folder, f'n_fourier{n_fourier}', f'T_anal_abb_x{xSIZE}y{ySIZE}.npy')
    else:
        dump_fname = os.path.join(upper_folder, f'n_fourier{n_fourier}', f'T_anal_eq_x{xSIZE}y{ySIZE}.npy')

    if os.path.isfile(dump_fname):
        print(f'{dump_fname} found, loading results from disc')
        (xx, yy, T_anal) = np.load(dump_fname)
        # T_anal = np.load(dump_fname)
    else:
        print(f'{dump_fname} not found, starting calculations')
        xx, yy, T_anal = analytical_laplace_2d(anal_input)
        fourier_folder = os.path.dirname(dump_fname)
        if not os.path.exists(fourier_folder):
            os.makedirs(fourier_folder)
        np.save(dump_fname, (xx, yy, T_anal))

    return xx, yy, T_anal


def make_anal_plot(xx, yy, zz):
    print("---------- PLOTTING -------------")
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xx, yy, zz, cmap=colormap.coolwarm,
                           linewidth=1, antialiased=False)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_zlim(0., 1.05)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
