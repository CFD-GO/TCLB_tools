import numpy as np
import time
from matplotlib import pyplot    #and the useful plotting library
import matplotlib.pylab as pylab
import os

from Benchmarks.Zaraza.SIR_1D_model import SIR_1D_FD_Peng, WSIR_1D_FD, SIR_1D_convolve_FD, SIR_1D_test
from Benchmarks.Zaraza.SIR_1D_model import make_wsir_plot_1D, compare_sir_vs_wsir_plot, make_hist_of_diffusivity_plot


nx = 129
domain_length = 64
dx = domain_length / (nx-1)
xspace = np.linspace(0, domain_length, nx)

r0 = 5.5  # infectious radius
beta_sir = 3.01  # the average number of contacts per person per time
gamma_sir = 1/2.8  # 1 over days to recovery

beta_W = 1e3

total_time = 1e0 # OK
dt = 1e-5 # OK

# total_time = 1e-0 # To crash
# dt = 1.25e-2 # To crash


ntimesteps = int(total_time / dt)

I_IC = np.ones(nx)*0.01                # numpy function ones()
I_IC[int((nx-1)/4):int(nx/2 + 1)] = 0.05  # setting u = 2 between 0.5 and 1 as per our I.C.s
S_IC = np.ones(nx) - I_IC
R_IC = np.zeros(nx)

N = S_IC + I_IC + R_IC
# make_wsir_plot_1D(S_IC, I_IC, R_IC, xspace, 0, 0, 'SIR IC')
start = time.process_time()

S, I, R, hist_of_diffusivity = SIR_1D_FD_Peng(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, ntimesteps, dt)
make_hist_of_diffusivity_plot(hist_of_diffusivity, ntimesteps, dt, title='SIR hist of diffusivity')
make_wsir_plot_1D(S, I, R, xspace, ntimesteps, dt, 'SIR 1D')

# Sw, Iw, Rw, Ww = WSIR_1D_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, ntimesteps, dt, beta_W)
# make_wsir_plot_1D(Sw, Iw, Rw, xspace, ntimesteps, dt, 'WSIR 1D', Ww)

# S, I, R = SIR_1D_convolve_FD(S_IC, I_IC, R_IC, nx, dx, r0, beta_sir, gamma_sir, ntimesteps, dt)
# make_wsir_plot_1D(S, I, R, xspace, ntimesteps, dt, 'SIR convolve 1D')

# compare_sir_vs_wsir_plot((S, I, R), (Sw, Iw, Rw, Ww), beta_W, xspace, ntimesteps, dt)


print('\n\n Done in %s [s].'
      % str(time.process_time() - start))
