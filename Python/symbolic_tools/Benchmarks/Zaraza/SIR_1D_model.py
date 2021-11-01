import numpy as np               #loading our favorite library
from numba import jit
from matplotlib import pyplot    #and the useful plotting library
import matplotlib.pylab as pylab
import os
from numba import vectorize, float64


def make_hist_of_diffusivity_plot(hist_of_diffusivity, nt, dt, title):
    how_many_slices = 100
    # sub_hist = hist_of_diffusivity[::int(nt / how_many_slices), :]  # take slice at each 100 timestep
    sub_hist = hist_of_diffusivity[:100, :]
    # sub_hist = hist_of_diffusivity[-100:, :]
    # sub_hist = hist_of_diffusivity
    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
             'axes.labelsize': 'xx-large',
             'axes.titlesize':'xx-large',
             'xtick.labelsize':'xx-large',
             'ytick.labelsize':'xx-large'}
    pylab.rcParams.update(params)
    axes = pyplot.gca()
    pos = pyplot.imshow(sub_hist, cmap='coolwarm',  interpolation='none')

    pyplot.title(f'{title} @ nt: {nt} dt {dt}')
    pyplot.xlabel(r'$x$')
    pyplot.ylabel(f'Time x {how_many_slices}')

    fig_name = f'plots/{title}.png'
    fig = pyplot.gcf()  # get current figure
    fig.colorbar(pos, ax=axes)
    pyplot.show()
    fig.savefig(fig_name, bbox_inches='tight')


def make_wsir_plot_1D(s, i, r, x, nt, dt, title, w=None):
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

    pyplot.xlabel(r'$x$')
    pyplot.ylabel(r'No of people')

    if w is not None:
        pyplot.plot(x, w, color="blue", linewidth=2, linestyle=":", label='W')

    axes.set_ylim([-0.05, 1.05])

    pyplot.legend()
    pyplot.title(f'{title} @ nt: {nt} dt {dt}')

    if not os.path.exists('plots'):
        os.makedirs('plots')

    fig_name = f'plots/{title}.png'
    fig = pyplot.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    pyplot.show()


def compare_sir_vs_wsir_plot(SIR, WSIR, beta_W, xspace, ntimesteps, dt):
    S, I, R = SIR
    Sw, Iw, Rw, Ww = WSIR

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (16, 10),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    axes = pyplot.gca()
    pyplot.plot(xspace, S, color="green", linewidth=2, marker='', markevery=2, label='S')
    pyplot.plot(xspace, I, color="red", linewidth=2, label='I')
    pyplot.plot(xspace, R, color="purple", linewidth=2, label='R')
    pyplot.plot(xspace, S + I + R, color="black", linewidth=2, label='Total population')

    pyplot.plot(xspace, Sw, color="green", linewidth=2, linestyle="", marker='o', markersize=4, markevery=3, label='S*')
    pyplot.plot(xspace, Iw, color="red", linewidth=2, linestyle="", marker='o', markersize=4, markevery=3, label='I*')
    pyplot.plot(xspace, Rw, color="purple", linewidth=2, linestyle="", marker='o', markersize=4, markevery=3,
                label='R*')
    pyplot.plot(xspace, Ww, color="blue", linewidth=2, linestyle="", marker='o', markersize=4, markevery=3, label='W*')

    pyplot.xlabel(r'$x$')
    pyplot.ylabel(r'No of people')

    axes.set_ylim([-0.05, 1.05])

    pyplot.legend()
    pyplot.title(f'SIR vs WSIR \n'
                 f'@ ntimesteps={ntimesteps:.2e}     dt={dt:.2e}     beta_W={beta_W:.2}')

    if not os.path.exists('plots'):
        os.makedirs('plots')
    fig_name = f'plots/SIRvsWSIR_betaW{beta_W:.2e}.png'
    fig = pyplot.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    pyplot.show()


@jit(cache=True, nopython=True)
def SIR_1D_FD_Peng(S, I, R, nx, dx, r0, beta_sir, gamma_sir, nt, dt):
    N = S + I + R

    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    hist_of_diffusivity = np.zeros((nt, nx), dtype=np.float64)
    for n in range(nt):  # iterate through time
        lap_I = (I[l_ind] - 2 * I[c_ind] + I[r_ind]) / dx ** 2
        qS2I_spatial = (r0 * r0 / 8.) * lap_I
        # qS2I_spatial = np.zeros(nx)

        hist_of_diffusivity[n] = beta_sir * S * qS2I_spatial
        # hist_of_diffusivity[n] = lap_I
        qS2I = dt * beta_sir * S * (qS2I_spatial + I) / N
        qI2R = dt * gamma_sir * I
        S = S - qS2I
        I = I + qS2I - qI2R
        R = R + qI2R



        # lap_S = (S[l_ind] - 2 * S[c_ind] + S[r_ind]) / dx ** 2
        # qS2S_spatial = (r0 * r0 / 8.) * lap_S
        # # qS2S_spatial = np.zeros(nx)
        # S = S - dt * beta_sir * (qS2S_spatial ) / N
        # I = I + dt * beta_sir * (qS2I_spatial ) / N #- qI2R
        # R = R #+ qI2R

        # if n % 10000 == 0:
        #     # Courant condition is: dt =< dx*dx/(2 * nu)
        #     print(f"calc_Courant_min = {min(C*dt/dx**2)} \t"
        #           f"dt_min to satisfy Courant condition= {dx**2/max(2*C)} \t"
        #           f"Growth factor S2I = {max(1 + dt * beta_sir * S * I / N)} \t"
        #           f"Growth factor I2R = {max(1 + dt * gamma_sir * I)}."   # Euler is stable if < 1
        #           )

    return S, I, R, hist_of_diffusivity

@jit(cache=True, nopython=True)
def SIR_1D_test(S, I, R, nx, dx, r0, beta_sir, gamma_sir, nt, dt):
    N = S + I + R

    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    hist_of_diffusivity = np.zeros((nt, nx), dtype=np.float64)
    for n in range(nt):  # iterate through time
        lap_I = (I[l_ind] - 2 * I[c_ind] + I[r_ind]) / dx ** 2
        hist_of_diffusivity[n] = lap_I
        lap_S = (S[l_ind] - 2 * S[c_ind] + S[r_ind]) / dx ** 2

        S = S + dt * beta_sir * (lap_S )
        I = I + dt * beta_sir * (lap_I )
        R = R

        # if n % 10000 == 0:
        #     # Courant condition is: dt =< dx*dx/(2 * nu)
        #     print(f"calc_Courant_min = {min(C*dt/dx**2)} \t"
        #           f"dt_min to satisfy Courant condition= {dx**2/max(2*C)} \t"
        #           f"Growth factor S2I = {max(1 + dt * beta_sir * S * I / N)} \t"
        #           f"Growth factor I2R = {max(1 + dt * gamma_sir * I)}."   # Euler is stable if < 1
        #           )

    return S, I, R, hist_of_diffusivity

# @jit(cache=True, nopython=True)
def SIR_1D_convolve_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, nt, dt):
    S = S_IC.copy()  # our placeholder array, to advance the solution in time
    I = I_IC.copy()
    R = R_IC.copy()
    N = S + I + R

    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    domain_length = dx*(nx-1)
    xspace = np.linspace(0, domain_length, nx)

    def get_gaussian(x, alfa, t):
        g = -(x-domain_length/2.)**2
        g /=(4*alfa*t)
        g = np.exp(g)
        g /= np.sqrt(4*np.pi*alfa*t)
        g *= domain_length/(nx-1) # normalize --> sum(g)=1
        return g

    window = get_gaussian(xspace, alfa=r0/2, t=1)
    # window[64] = 1

    # window = np.zeros(nx)
    # window[65] =1./3
    # window[64] = 1. / 3
    # window[63] = 1. / 3

    for n in range(nt):  # iterate through time
        # lap_I = (I[l_ind] - 2 * I[c_ind] + I[r_ind]) / dx ** 2
        # qS2I_spatial = (r0 * r0 / 8.) * lap_I
        # qS2I_spatial = np.zeros(nx)

        # qS2I = dt * beta_sir * S * (qS2I_spatial + I) / N
        qS2I = dt * beta_sir * S * (np.convolve(window, I, 'same')) / N
        # qS2I = dt * beta_sir * S * conv_circ(I, window,) / N
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
def WSIR_1D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, nt, dt, beta_W=1e2):
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

        qW = dt * beta_W * (qW_spatial + (I - W))
        qS2I = dt * beta_sir * S * W/N
        qI2R = dt * gamma_sir * I

        W = W + qW
        S = S - qS2I
        I = I + qS2I - qI2R
        R = R + qI2R

        N = S + I + R  # TODO: uzywać przeliczonego N czy N z IC?

    return S, I, R, W

