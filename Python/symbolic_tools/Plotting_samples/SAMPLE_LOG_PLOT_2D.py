from matplotlib import pyplot as plt
import matplotlib.ticker
import numpy as np
import os

coeff = np.array([1, 2, 4, 8, 16, 32])
x = 1*coeff
initial_error = 1
y_1st = initial_error/coeff
y_2nd = initial_error/(coeff*coeff)

# ------------------------------------ PLOT ------------------------------------
if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/sample_loglog_plot.png'

fig1, ax1 = plt.subplots(figsize=(14, 8))
plt.rcParams.update({'font.size': 14})
ax1.plot(x, y_1st,
         color="black", marker=">", markevery=1, markersize=8, linestyle=":", linewidth=2,
         label=r'first order convergence, $\mathcal{O}(n)$')

ax1.plot(x, y_2nd,
         color="black", marker="D", markevery=1, markersize=8, linestyle="-", linewidth=2,
         label=r'second order convergence, $\mathcal{O}(n^2)$')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xticks(x)

ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.title(f'Sample log plot\n '
          r'$x_{range}$' + f'={(x.max()-x.min())}'
          f'; \t'
          r'$x_{step}$' + f'={len(x):.4f}')
plt.xlabel(r'$x_{label}$', fontsize=18)
plt.ylabel(r'$y_{label}$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=12)
plt.legend()
plt.grid(True,  which="both")
fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()
# plt.close(fig)  # close the figure
