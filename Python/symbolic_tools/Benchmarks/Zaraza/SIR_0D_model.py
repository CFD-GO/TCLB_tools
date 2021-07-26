# inspirations:
# http://gabgoh.github.io/COVID/index.html?fbclid=IwAR2VJC06vWXCwTKuwC0-DdkhqJGjzlIdZAWFtsII5wn7VWbItlTXONmChNQ
# https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy.integrate import solve_ivp


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


# CONSTANTS
beta = 0.71  # number of contacts per day
gamma = 1/3.2  # 1 over days to recovery

N = 1e6  # Size of population [no of people].

# INITIAL CONDItIONS
initial_susceptible = 0.99*N  # initial number of susceptible individuals in population.
initial_infections = 0.01*N  # initial number of infected individuals in population.
initial_removed = 0  # initial number of removed (recovered) individuals in population.
IC = np.array([initial_susceptible, initial_infections, initial_removed])

days_to_simulate = 30
sol = solve_ivp(SIR_0D,
                [0, days_to_simulate],
                IC,
                method='RK45',
                args=[beta, gamma, N],
                dense_output=True)

t = np.linspace(0, days_to_simulate, 1000)
z = sol.sol(t)

S, I, R = z

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (14, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
axes = plt.gca()
plt.plot(t, S,
         color="green", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='Susceptible')
plt.plot(t, I,
         color="red", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='Infected')
plt.plot(t, R,
         color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='Removed')


plt.xlabel('t [days]')
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.ylabel('No of people')
plt.title('SIR Epidemic Calculator')
plt.legend()
plt.grid()
plt.show()

