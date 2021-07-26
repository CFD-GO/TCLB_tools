import numpy as np               #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library

nx = 128
domain_length = 64
dx = domain_length / (nx-1)
xspace = np.linspace(0, domain_length, nx)

r0 = 0.025  # infectious radius
beta_sir = 3.01  # the average number of contacts per person per time
gamma_sir = 1/3.2  # 1 over days to recovery

total_time = 1

I_IC = np.ones(nx)*0.05                  # numpy function ones()
I_IC[int((nx-1)/4):int(nx/2 + 1)] = 0.1  # setting u = 2 between 0.5 and 1 as per our I.C.s
S_IC = np.ones(nx) - I_IC
R_IC = np.zeros(nx)

N = S_IC + I_IC + R_IC

#
# pyplot.plot(xspace, S_IC, color="green", label='Susceptible', marker='x', markevery=2)
# pyplot.plot(xspace, I_IC, color="red", label='Infected')
# pyplot.plot(xspace, R_IC, color="purple", label='Recovered')
# pyplot.plot(xspace, N, color="black", label='Total population')
#
# pyplot.legend()
# pyplot.title('Initial Condition')
# pyplot.show()


def calc_diffusion_FD(S_IC, I_IC, R_IC, nx, r0, beta_sir, gamma_sir, total_time):
    S = S_IC.copy()
    # Sn = S_IC.copy()  # our placeholder array, un, to advance the solution in time

    I = I_IC.copy()
    # In = I_IC.copy()  # our placeholder array, un, to advance the solution in time

    R = R_IC.copy()
    # Rn = R_IC.copy()  # our placeholder array, un, to advance the solution in time

    C0 = beta_sir * r0 * r0 / 8.
    Courant = 0.2  # the condition is: dt =< dx*dx/(2 * nu)
    #     dt = sigma * dx**2 / nu

    dt = 1e-6
    # nt = int(total_time/dt)
    nt = 10000
    # qSI_spatial = np.zeros(nx)  # placeholder for coupling terms
    # qSI = np.zeros(nx)

    c_ind = np.arange(0, nx)
    l_ind = np.roll(c_ind, -1)
    r_ind = np.roll(c_ind, 1)

    for n in range(nt):  # iterate through time
        In = I.copy()  # copy the existing values of I into In
        Sn = S.copy()
        # # Rn = R.copy()

        qSI = dt * beta_sir * Sn * In / N

        C = C0 * Sn * dt / N
        lap_I = (In[l_ind] - 2 * In[c_ind] + In[r_ind]) / dx ** 2
        qSI_spatial = C * lap_I

        S = S - qSI - qSI_spatial
        I = I + qSI + qSI_spatial
        R = R + dt * gamma_sir * I

        # for i in range(0, nx - 1): # the slow way - (and bad BC)
        #     qSI = dt * beta_sir * Sn[i] * In[i] / N[i]
        #     # qSI[i] = 0
        #
        #     lap_I = (In[i - 1] - 2*In[i] + In[i + 1])/dx**2
        #     C = C0 * Sn[i] * dt /N[i]
        #     qSI_spatial = C * lap_I
        #     # qSI_spatial[i] = 0
        #
        #     S[i] = S[i] - qSI - qSI_spatial
        #     I[i] = I[i] + qSI + qSI_spatial
        #     R[i] = R[i] + dt * gamma_sir * I[i]

    return S, I, R


S, I, R = calc_diffusion_FD(S_IC, I_IC, R_IC, nx, r0, beta_sir, gamma_sir, total_time)

pyplot.plot(xspace, S, color="green", label='Susceptible', marker='x', markevery=2)
pyplot.plot(xspace, I, color="red", label='Infected')
pyplot.plot(xspace, R, color="purple", label='Recovered')
pyplot.plot(xspace, S + I + R, color="black", label='Total population')

pyplot.legend()
pyplot.title('Initial Condition')
pyplot.show()
