# References:
# https://www.wire.tu-bs.de/lehre/ws15/pde1/lecture_2.pdf
# https://www.wire.tu-bs.de/lehre/ws15/pde1/lecture_3.pdf


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
import sympy as sp
from sympy.abc import x, y, t
from sympy import pretty_print
from sympy import fourier_series, pi, cos, sin, exp, pi, integrate
from typing import Callable
# from dataclasses import dataclass
from joblib import Parallel, delayed
import multiprocessing
import time
from typing import NamedTuple


# def analytical_laplace_2D(input):
def analytical_laplace_2D():
    # xySIZE = 1
    # step = 0.1
    # # my_fun = -4 * x / (x - xySIZE)/(xySIZE*xySIZE)
    # my_fun = x / (1 - x)

    x_high = 1
    y_high = 1
    step = 0.1
    my_fun = x / (1 - x)
    x_low = 0
    y_low = 0
    n_fourier_terms = 25

    lim = (x, x_low, x_high)
    Lx = x_high - x_low
    Ly = y_high - y_low
    u_sol = 0
    c = [0]

    start = time.process_time()
    print("---------- Calculating Fourier coeff -------------")
    for k in range(1, n_fourier_terms):  # skip zero
        result = 2 / (Lx * sp.sinh(k * pi * Ly / Lx)) * integrate(my_fun * sin(k * pi * x / Lx), lim)
        c.append(result)
        u_sol += c[k] * sp.sinh(k * pi * y / Lx) * sp.sin(k * pi * x / Lx)

    print("---------- Calculating Values -------------")

    nx = int((x_high - x_low) / step)
    ny = int((y_high - y_low) / step)
    x_grid = np.linspace(x_low, x_high, nx)
    y_grid = np.linspace(y_low, y_high, ny)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = np.zeros((ny, nx))

    # lim = (x, x_low, x_high)
    # Lx = x_high - x_low
    # Ly = y_high - y_low
    # u_sol = 0
    # c = [0]

    # start = time.process_time()
    # print("---------- Calculating Fourier coeff -------------")
    # for k in range(1, input.n_fourier_terms):  # skip zero
    #     result = 2 / (Lx * sp.sinh(k * pi * Ly / Lx)) * integrate(input.my_fun * sin(k * pi * x / Lx), lim)
    #     c.append(result)
    #     u_sol += c[k] * sp.sinh(k * pi * y / Lx) * sp.sin(k * pi * x / Lx)
    #
    # print("---------- Calculating Values -------------")
    #
    # nx = int((input.x_high - input.x_low) / input.step)
    # ny = int((input.y_high - input.y_low) / input.step)
    # x_grid = np.linspace(input.x_low, input.x_high, nx)
    # y_grid = np.linspace(input.y_low, input.y_high, ny)
    # X, Y = np.meshgrid(x_grid, y_grid)
    # Z = np.zeros((ny, nx))

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
    for i in range(ny):
        # if i % ny/10 == 0:
        print(f"=== Done i/ny: {i}/{ny}  ===")

        for j in range(nx):
            print(f"Doing i/ny: {i}/{ny} \t j/nx: {j}/{nx}")
            # Z[i][j] = my_fun.subs({'x': X[i][j]})
            Z[i][j] = u_sol.subs({'x': X[i][j], 'y': Y[i][j]})

    print(f'\n\n Done in {time.process_time() - start} [s].')
    return X, Y, Z


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
    # ax.set_zlim(-1.01, 1.01)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


@dataclass
class InputForLaplace2DAnalytical:
    x_high: int
    y_high: int
    step: float
    my_fun: Callable
    x_low: int = 0
    y_low: int = 0
    n_fourier_terms: int = 5
