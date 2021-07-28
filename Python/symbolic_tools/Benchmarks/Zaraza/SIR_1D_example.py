import numpy as np
import time
from matplotlib import pyplot    #and the useful plotting library
import matplotlib.pylab as pylab
import os
from Benchmarks.Zaraza.SIR_1D_model import SIR_1D_FD, WSIR_1D_FD
from Benchmarks.Zaraza.SIR_1D_model import plotuj_wsira_1D



nx = 128
domain_length = 64
dx = domain_length / (nx-1)
xspace = np.linspace(0, domain_length, nx)

r0 = 15.5  # infectious radius
beta_sir = 3.01  # the average number of contacts per person per time
gamma_sir = 1/3.2  # 1 over days to recovery

beta_W = 1e4

total_time = 1e0
dt = 1e-7
ntimesteps = int(total_time / dt)

I_IC = np.ones(nx)*0.01                # numpy function ones()
I_IC[int((nx-1)/4):int(nx/2 + 1)] = 0.05  # setting u = 2 between 0.5 and 1 as per our I.C.s
S_IC = np.ones(nx) - I_IC
R_IC = np.zeros(nx)

N = S_IC + I_IC + R_IC
# plotuj_wsira(S_IC, I_IC, R_IC, xspace, time_spot=0)
start = time.process_time()

S, I, R = SIR_1D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, ntimesteps, dt)
# plotuj_wsira_1D(S, I, R, xspace, nt, dt, 'SIR 1D')

Sw, Iw, Rw, Ww = WSIR_1D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, ntimesteps, dt, beta_W)
# plotuj_wsira_1D(Sw, Iw, Rw, xspace, nt, dt, 'WSIR 1D', Ww)

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
pyplot.plot(xspace, Iw, color="red", linewidth=2, linestyle="",  marker='o', markersize=4, markevery=3, label='I*')
pyplot.plot(xspace, Rw, color="purple", linewidth=2, linestyle="",  marker='o', markersize=4, markevery=3, label='R*')
pyplot.plot(xspace, Ww, color="blue", linewidth=2, linestyle="",  marker='o', markersize=4, markevery=3, label='W*')

axes.set_ylim([-0.05, 1.05])

pyplot.legend()
pyplot.title(f'SIR vs WSIR \n'
             f'@ ntimesteps={ntimesteps:.2e}     dt={dt:.2e}     beta_W={beta_W:.2}')


if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/SIRvsWSIR_beta{beta_W:.2e}.png'
fig = pyplot.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
pyplot.show()

print('\n\n Done in %s [s].'
      % str(time.process_time() - start))
