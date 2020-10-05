# inspirations:
# http://gabgoh.github.io/COVID/index.html?fbclid=IwAR2VJC06vWXCwTKuwC0-DdkhqJGjzlIdZAWFtsII5wn7VWbItlTXONmChNQ
# https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy.integrate import solve_ivp


def SEIR(t, z, R_0, T_inc, T_inf, N):
    """
     # Susceptible → Exposed → Infected → Removed
    :param t: time [days]
    :param z: Susceptible, Exposed, Infected, Removed
    :param T_inf: Duration patient is infectious
    :param R_0: Basic Reproduction Number. Measure of contagiousness: the number of secondary infections each infected individual produces.
    :param T_inc: Length of incubation period
    :return: derivatives [dS, dE, dI, dR]
    """
    beta = R_0 /T_inf

    S, E, I, R = z
    dS = -beta*I*S/N
    dE = beta*I*S/N-E/T_inc
    dI = E/T_inc - I/T_inf
    dR = I/T_inf

    # if round(t*1, 0) % 1 == 0:
    #     print(f"solving: {t} [day]")

    return [dS, dE, dI, dR]


# CONSTANTS
R0 = 2.2  # Basic Reproduction Number - the number of secondary infections each infected individual produces.
T_inc = 5.2  # Length of incubation period [days]
T_inf = 2.9  # Duration patient is infectious [days]
N = 1e6  # Size of population [no of people].

# INITIAL CONDItIONS
initial_susceptible = 0.999 * N  # initial number of susceptible individuals in population.
initial_exposed = 0.001*N  # initial number of exposed individuals in population.
initial_infections = 0*N  # initial number of infected individuals in population.
initial_removed = 0  # initial number of removed (recovered) individuals in population.
IC = np.array([initial_susceptible, initial_exposed, initial_infections, initial_removed])

days_to_simulate = 150
sol = solve_ivp(SEIR,
                [0, days_to_simulate],
                IC,
                method='RK45',
                args=[R0, T_inc, T_inf, N],
                dense_output=True)

t = np.linspace(0, days_to_simulate, 1000)
z = sol.sol(t)

S, E, I, R = z

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
plt.plot(t, E,
         color="blue", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='Exposed')
plt.plot(t, I,
         color="red", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='Infected')
plt.plot(t, R,
         color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label='Removed')

plt.xlabel('t [days]')
plt.ylabel('No of people')
plt.title('SEIR Epidemic Calculator')
plt.legend()
plt.grid()
plt.show()
