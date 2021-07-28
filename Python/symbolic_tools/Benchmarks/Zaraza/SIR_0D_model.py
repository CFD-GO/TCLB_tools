# inspirations:
# http://gabgoh.github.io/COVID/index.html?fbclid=IwAR2VJC06vWXCwTKuwC0-DdkhqJGjzlIdZAWFtsII5wn7VWbItlTXONmChNQ
# https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab


import numpy as np
def plotuj_wsira_0D(s, i, r, t, title, w=None):
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


def WSIR_0D(t, z, beta, gamma, beta_LLW, N):
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

    dWdt = beta_LLW * (I-W)
    dSdt = -beta*W*S/N
    dIdt = beta*W*S/N - I*gamma
    dRdt = I*gamma
    return [dSdt, dIdt, dRdt, dWdt]


