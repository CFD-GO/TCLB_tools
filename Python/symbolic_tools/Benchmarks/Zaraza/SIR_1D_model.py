import numpy as np               #loading our favorite library
from numba import jit
from matplotlib import pyplot    #and the useful plotting library
import matplotlib.pylab as pylab

from numba import vectorize, float64
# from Benchmarks.Zaraza.SIR_0D_model import SIR_0D

def plotuj_wsira_1D(s, i, r, x, nt, dt, title, w=None ):
    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
             'axes.labelsize': 'xx-large',
             'axes.titlesize':'xx-large',
             'xtick.labelsize':'xx-large',
             'ytick.labelsize':'xx-large'}
    pylab.rcParams.update(params)
    axes = pyplot.gca()
    pyplot.plot(x, s, color="green", linewidth=2, marker='', markevery=2, label='Susceptible')
    pyplot.plot(x, i, color="red", linewidth=2, label='Infected')
    pyplot.plot(x, r, color="purple", linewidth=2, label='Recovered')
    pyplot.plot(x, s + i + r, color="black", linewidth=2, label='Total population')

    if w is not None:
        pyplot.plot(x, w, color="blue", linewidth=2, linestyle=":", label='W')

    axes.set_ylim([-0.05, 1.05])

    pyplot.legend()
    pyplot.title(f'{title} @ nt: {nt} dt {dt}')
    pyplot.show()


@jit(cache=True, nopython=True)
def SIR_1D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, nt, dt):
    S = S_IC.copy()  # our placeholder array, to advance the solution in time
    I = I_IC.copy()
    R = R_IC.copy()
    N = S + I + R

    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    for n in range(nt):  # iterate through time
        lap_I = (I[l_ind] - 2 * I[c_ind] + I[r_ind]) / dx ** 2
        qS2I_spatial = (r0 * r0 / 8.) * lap_I
        # qS2I_spatial = np.zeros(nx)

        qS2I = dt * beta_sir * S * (qS2I_spatial + I) / N
        qI2R = dt * gamma_sir * I
        S = S - qS2I
        I = I + qS2I - qI2R
        R = R + qI2R

        N = S + I + R
        # if n % 10000 == 0:
        #     # Courant condition is: dt =< dx*dx/(2 * nu)
        #     print(f"calc_Courant_min = {min(C*dt/dx**2)} \t"
        #           f"dt_min to satisfy Courant condition= {dx**2/max(2*C)} \t"
        #           f"Growth factor S2I = {max(1 + dt * beta_sir * S * I / N)} \t"
        #           f"Growth factor I2R = {max(1 + dt * gamma_sir * I)}."   # Euler is stable if < 1
        #           )

    return S, I, R



@jit(cache=True, nopython=True)
def WSIR_1D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, nt, dt, beta_LLW=1e2):
    S = S_IC.copy()  # our placeholder array, to advance the solution in time
    I = I_IC.copy()
    R = R_IC.copy()
    W = np.zeros(nx)

    N = S + I + R
    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    for n in range(nt):  # iterate through time
        lap_W = (W[l_ind] - 2 * W[c_ind] + W[r_ind]) / dx ** 2
        qW_spatial = (r0 * r0 / 8.)*lap_W
        # qW_spatial = np.zeros(nx)

        qW = dt * beta_LLW * (qW_spatial + (I - W))
        qS2I = dt * beta_sir * S * W/N
        qI2R = dt * gamma_sir * I

        W = W + qW
        S = S - qS2I
        I = I + qS2I - qI2R
        R = R + qI2R

        N = S + I + R  # TODO: uzywać przeliczonego N czy N z IC?

    return S, I, R, W

