# inspirations:
# http://gabgoh.github.io/COVID/index.html?fbclid=IwAR2VJC06vWXCwTKuwC0-DdkhqJGjzlIdZAWFtsII5wn7VWbItlTXONmChNQ
# https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from numba import jit

import numpy as np
def make_wsir_plot_0D(s, i, r, t, title, w=None):
    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)
    axes = plt.gca()
    plt.plot(t, s,
             color="green", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
             label='Susceptible')
    plt.plot(t, i,
             color="red", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
             label='Infected')
    plt.plot(t, r,
             color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
             label='Removed')
    if w is not None:
        plt.plot(t, w,
                 color="blue", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
                 label='W')

    plt.xlabel('t [days]')
    # plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.ylabel('No of people')
    plt.title(f'{title}')
    plt.legend()
    plt.grid()
    plt.show()


def SIR_0D(t, z, beta, gamma, N):
    """
     # Susceptible → Infected → Removed
    :param t: time [days]
    :param z: Susceptible, Exposed, Infected, Removed
    :param beta: average number of contacts per day for each infected individual
    :param gamma: Between I and R, the transition rate is γ
    (simply the frequency of recoveries, that is, number of recovered or dead during one day
    divided by the total number of infected on that same day, supposing "day" is the time unit).
    If the duration of the infection is denoted D, then γ = 1/D.
    :return: derivatives [dS, dI, dR]
    """

    S, I, R = z
    dSdt = -beta*I*S/N
    dIdt = beta*I*S/N - I*gamma
    dRdt = I*gamma
    return [dSdt, dIdt, dRdt]


def WSIR_0D(t, z, beta, gamma, beta_W, N):
    """
     # Susceptible → Infected → Removed
    :param t: time [days]
    :param z: Susceptible, Exposed, Infected, Removed
    :param beta: average number of contacts per day for each infected individual
    :param gamma: Between I and R, the transition rate is γ
    (simply the frequency of recoveries, that is, number of recovered or dead during one day
    divided by the total number of infected on that same day, supposing "day" is the time unit).
    If the duration of the infection is denoted D, then γ = 1/D.
    :return: derivatives [dS, dI, dR]
    """

    S, I, R, W = z

    dWdt = beta_W * (I - W)
    dSdt = -beta*W*S/N
    dIdt = beta*W*S/N - I*gamma
    dRdt = I*gamma
    return [dSdt, dIdt, dRdt, dWdt]

# solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, dense_output=False,
#               events=None, vectorized=False, args=None, **options):
@jit(cache=True, nopython=True)
def SIR_0D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, nt, dt):
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
        qS2I_spatial = np.zeros(nx)

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