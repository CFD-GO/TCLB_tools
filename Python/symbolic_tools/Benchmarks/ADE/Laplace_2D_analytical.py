
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


@dataclass
class InputForLaplace2DAnalytical:
    x_high: int
    y_high: int
    step: float
    my_fun: Callable
    n_fourier_terms: int
    x_low: int = 0
    y_low: int = 0


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

    nx = int((input_config.x_high - input_config.x_low) / input_config.step)
    ny = int((input_config.y_high - input_config.y_low) / input_config.step)
    x_grid = np.linspace(input_config.x_low, input_config.x_high, nx, endpoint=False) + 0.5
    y_grid = np.linspace(input_config.y_low, input_config.y_high, ny, endpoint=False) + 0.5
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
        print(f"=== Doing i/ny: {i}/{ny}  ===")
        for j in range(nx):
            # print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
            # zz[i][j] = my_hacked_fun.subs({'x': xx[i][j]})
            zz[i][j] = u_sol.subs({'x': xx[i][j], 'y': yy[i][j]})

    print(f'\n\n Done in {time.process_time() - start} [s].')
    return xx, yy, zz


def make_anal_plot(xx, yy, zz):
    print("---------- PLOTTING -------------")

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Plot the surface.
    surf = ax.plot_surface(xx, yy, zz, cmap=colormap.coolwarm,
                           linewidth=1, antialiased=False)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Customize the z axis.
    ax.set_zlim(0., 1.05)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
