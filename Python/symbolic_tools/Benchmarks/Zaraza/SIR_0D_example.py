import numpy as np

from scipy.integrate import solve_ivp

from Benchmarks.Zaraza.SIR_0D_model import SIR_0D, WSIR_0D
from Benchmarks.Zaraza.SIR_0D_model import plotuj_wsira_0D


# CONSTANTS
beta = 0.71  # number of contacts per day
gamma = 1/3.2  # 1 over days to recovery

N = 1e6  # Size of population [no of people].

# INITIAL CONDItIONS
initial_susceptible = 0.99*N  # initial number of susceptible individuals in population.
initial_infections = 0.01*N  # initial number of infected individuals in population.
initial_removed = 0  # initial number of removed (recovered) individuals in population.

initial_W = 0
beta_LLW = 1e2

IC = np.array([initial_susceptible, initial_infections, initial_removed])

days_to_simulate = 30
t = np.linspace(0, days_to_simulate, 1000)


sol = solve_ivp(SIR_0D,
                [0, days_to_simulate],
                IC,
                method='RK45',
                args=[beta, gamma, N],
                dense_output=True)
z = sol.sol(t)
S, I, R = z

plotuj_wsira_0D(S, I, R, t, 'SIR Epidemic Calculator', w=None)


ICw = np.array([initial_susceptible, initial_infections, initial_removed, initial_W])
sol = solve_ivp(WSIR_0D,
                [0, days_to_simulate],
                ICw,
                method='RK45',
                args=[beta, gamma, beta_LLW, N],
                dense_output=True)
zw = sol.sol(t)
Sw, Iw, Rw, Ww = zw
plotuj_wsira_0D(Sw, Iw, Rw, t, 'WSIR Epidemic Calculator', w=Ww)

print("Done.")
