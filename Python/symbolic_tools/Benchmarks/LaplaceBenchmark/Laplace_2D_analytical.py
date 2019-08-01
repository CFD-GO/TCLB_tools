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

import time

import os
import pwd


def peel_the_skin(some_2d_array):
    # clip the from each side, so that both EQ and ABB scheme would have the same 'measurement' nodes
    some_2d_array = np.delete(some_2d_array, 0, axis=0)
    some_2d_array = np.delete(some_2d_array, 0, axis=1)

    # n_rows, n_columns = some_2d_array.shape
    some_2d_array = np.delete(some_2d_array, (some_2d_array.shape[0] - 1), axis=0)  # delete last row
    some_2d_array = np.delete(some_2d_array, (some_2d_array.shape[1] - 1), axis=1)  # delete last column

    return some_2d_array


def prepare_anal_data_new(ySIZE, xSIZE, folder, shall_recalculate_results=False):
    # two functions are needed - one for ABB (with zeros in 1 and 63) and second for EQ (zeros in 0.5 and 63.5)
    def analytical_laplace_2d_new(is_bc_wet_node_type, n_fourier_terms):
        x1 = 0
        y1 = 0
        # is_bc_wet_node_type = True
        if is_bc_wet_node_type:
            x2 = xSIZE - 1
            y2 = ySIZE - 1
        else:
            x2 = xSIZE - 2
            y2 = ySIZE - 2

        xm = 0.5 * (x1 + x2)
        ym = 1
        a = ym / ((xm - x1) * (xm - x2))
        # my_fun = a * (x - x1) * (x - x2)
        my_fun = sp.sin(sp.pi*x/x2)

        # Prepare grid
        if is_bc_wet_node_type:
            x_grid = np.linspace(start=x1, stop=x2, num=x2+1, endpoint=True)
            y_grid = np.linspace(start=y1, stop=y2, num=y2+1, endpoint=True)
        else:
            # although the BC is effectively imposed in a mid-node location,
            # the (non-physical) values are also computed at the nodes
            # to make both wet-node/link-wise reference solutions of the same size.
            # x_grid = np.linspace(start=x1, stop=x2+1, num=x2+2, endpoint=True)
            # y_grid = np.linspace(start=y1, stop=y2+1, num=y2+2, endpoint=True)

            x_grid = np.linspace(start=x1-0.5, stop=x2+0.5, num=x2+2, endpoint=True)
            y_grid = np.linspace(start=y1-0.5, stop=y2+0.5, num=y2+2, endpoint=True)

            # x_grid = np.linspace(start=x1, stop=x2+1., num=x2+2, endpoint=True)
            # y_grid = np.linspace(start=y1, stop=y2+1., num=y2+2, endpoint=True)

        xx, yy = np.meshgrid(x_grid, y_grid)
        zz = np.zeros((ySIZE, xSIZE))

        # Calculate anal solution
        lim = (x, x1, x2)
        Lx = x2 - x1
        Ly = y2 - y1
        u_sol = 0
        c = [0]

        start = time.process_time()
        print("---------- Calculating Fourier coeff -------------")
        for k in range(1, n_fourier_terms):  # skip zero
            print(f"=== Calculating k/n_fourier_terms: {k}/{n_fourier_terms}  ===")
            result = 2 / (Lx * sp.sinh(k * pi * Ly / Lx)) * integrate(my_fun * sin(k * pi * x / Lx), lim)
            c.append(result)
            if is_bc_wet_node_type:
                u_sol += c[k] * sp.sinh(k * pi * y / Lx) * sp.sin(k * pi * x / Lx)
            else:
                u_sol += c[k] * sp.sinh(k * pi * (y-0.5) / Lx) * sp.sin(k * pi * x / Lx)

        print("---------- Calculating values on the grid nodes -------------")
        # Row-major order is also known as the C order, as the C programming language uses it.
        # New NumPy arrays are by default in row-major order.
        # my_hacked_fun = -4 * x * (x - input_config.x_high)/(input_config.x_high*input_config.x_high)
        for i in range(ySIZE):
            print(f"=== Calculating i/ny: {i}/{ySIZE}  ===")
            for j in range(xSIZE):
                # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
                zz[i][j] = u_sol.subs({'x': xx[i][j], 'y': yy[i][j]})

        # once the calculations are done,
        # we fix ideas (for postprocessing only) by shifting the grid so it would coincide with the one from TCLB
        # if is_bc_wet_node_type:
        #     xx += 0.5
        #     yy += 0.5
        # else:
        #     xx += 1
        #     yy += 1

        print(f'\n\n Done in {time.process_time() - start} [s].')
        return xx, yy, zz

    n_fourier = 2
    upper_folder = os.path.dirname(folder)

    if 'abb' in folder:
        is_bc_wet_node_type = False
        dump_fname = os.path.join(upper_folder, f'n_fourier{n_fourier}', f'T_anal_abb_x{xSIZE}y{ySIZE}.npy')
    else:
        is_bc_wet_node_type = True
        dump_fname = os.path.join(upper_folder, f'n_fourier{n_fourier}', f'T_anal_eq_x{xSIZE}y{ySIZE}.npy')

    if os.path.isfile(dump_fname) and not shall_recalculate_results:
        print(f'{dump_fname} found, loading results from disc')
        (xx, yy, T_anal) = np.load(dump_fname)
    else:
        print(f'recalculating results')
        xx, yy, T_anal = analytical_laplace_2d_new(is_bc_wet_node_type, n_fourier)
        fourier_folder = os.path.dirname(dump_fname)
        if not os.path.exists(fourier_folder):
            os.makedirs(fourier_folder)
        np.save(dump_fname, (xx, yy, T_anal))

    return xx, yy, T_anal
