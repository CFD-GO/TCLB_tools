import numpy as np               #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
import matplotlib.pylab as pylab
from numba import jit
import time
from numba import vectorize, float64
# from Benchmarks.Zaraza.SIR_0D_model import SIR_0D

start = time.process_time()

# @jit

nx = 128
domain_length = 64
dx = domain_length / (nx-1)
xspace = np.linspace(0, domain_length, nx)

r0 = 15.5  # infectious radius
beta_sir = 3.01  # the average number of contacts per person per time
gamma_sir = 1/3.2  # 1 over days to recovery

total_time = 1e0
dt = 1e-3

I_IC = np.ones(nx)*0.05                  # numpy function ones()
I_IC[int((nx-1)/4):int(nx/2 + 1)] = 0.6  # setting u = 2 between 0.5 and 1 as per our I.C.s
S_IC = np.ones(nx) - I_IC
R_IC = np.zeros(nx)

N = S_IC + I_IC + R_IC

#
# pyplot.plot(xspace, S_IC, color="green", label='Susceptible', marker='x', markevery=2)
# pyplot.plot(xspace, I_IC, color="red", label='Infected')
# pyplot.plot(xspace, R_IC, color="purple", label='Recovered')
# pyplot.plot(xspace, S_IC + I_IC + R_IC, color="black", label='Total population')
#
# pyplot.legend()
# pyplot.title('Initial Condition')
# pyplot.show()

@jit(cache=True, nopython=True)
def SIR_1D_FD(S_IC, I_IC, R_IC, nx, r0, beta_sir, gamma_sir, total_time, dt):
    S = S_IC.copy()  # our placeholder array, to advance the solution in time
    I = I_IC.copy()
    R = R_IC.copy()
    N = S + I + R

    C0 = beta_sir * r0 * r0 / 8.
    nt = int(total_time/dt)
    # nt = 1000

    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    for n in range(nt):  # iterate through time
        qS2I = dt * beta_sir * S * I / N

        C = C0 * S * dt / N
        lap_I = (I[l_ind] - 2 * I[c_ind] + I[r_ind]) / dx ** 2
        qS2I_spatial = C * lap_I
        # qS2I_spatial = np.zeros(nx)

        qI2R = dt * gamma_sir * I
        S = S - qS2I - qS2I_spatial
        I = I + qS2I + qS2I_spatial - qI2R
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

S, I, R = SIR_1D_FD(S_IC, I_IC, R_IC, nx, r0, beta_sir, gamma_sir, total_time, dt)

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (14, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
axes = pyplot.gca()
pyplot.plot(xspace, S, color="green", label='Susceptible', marker='x', markevery=2)
pyplot.plot(xspace, I, color="red", label='Infected')
pyplot.plot(xspace, R, color="purple", label='Recovered')
pyplot.plot(xspace, S + I + R, color="black", label='Total population')
axes.set_ylim([-0.05, 1.05])

pyplot.legend()
pyplot.title(f'SIR 1D @ time: {total_time}')
pyplot.show()

print('\n\n Done in %s [s].'
      % str(time.process_time() - start))

# from scipy.integrate import solve_ivp
#
# S_IC_0D = S_IC[int((nx-1)/4 + 2)]
# I_IC_0D = I_IC[int((nx-1)/4 + 2)]
# R_IC_0D = R_IC[int((nx-1)/4 + 2)]
# IC_0D = np.array([S_IC_0D, I_IC_0D, R_IC_0D])
# N_0D = sum(IC_0D)
#
# sol = solve_ivp(SIR_0D,
#                 [0, total_time],
#                 IC_0D,
#                 method='RK45',
#                 args=[beta_sir, gamma_sir, N_0D],
#                 dense_output=True)
#
# t = np.linspace(0, total_time, 1000)
# z = sol.sol(t)
#
# S_0D, I_0D, R_0D = z
#
# import matplotlib.pyplot as plt
# import matplotlib.pylab as pylab
# params = {'legend.fontsize': 'xx-large',
#           'figure.figsize': (14, 8),
#          'axes.labelsize': 'xx-large',
#          'axes.titlesize':'xx-large',
#          'xtick.labelsize':'xx-large',
#          'ytick.labelsize':'xx-large'}
# pylab.rcParams.update(params)
# axes = plt.gca()
# plt.plot(t, S_0D,
#          color="green", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
#          label='Susceptible')
# plt.plot(t, I_0D,
#          color="red", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
#          label='Infected')
# plt.plot(t, R_0D,
#          color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
#          label='Removed')
#
#
# plt.xlabel('t [days]')
# # plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
# plt.ylabel('No of people')
# plt.title('SIR Epidemic Calculator')
# plt.legend()
# plt.grid()
# plt.show()
