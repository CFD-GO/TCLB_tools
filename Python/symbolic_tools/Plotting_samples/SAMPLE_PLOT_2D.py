import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import os
# -------------------- prepare dummy data --------------------

x_range = 1
step = 0.01
x = np.arange(0., x_range, step)
y = -4 * x * (x - x_range) / (x_range * x_range)

if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/sample_plot2D_param_b={x_range}.png'

# -------------------- make dummy plot --------------------
# plt.figure(figsize=(14, 8))
# plt.rcParams.update({'font.size': 14})

# https://matplotlib.org/api/font_manager_api.html#matplotlib.font_manager.FontProperties.set_size
params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (14, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)



axes = plt.gca()
plt.plot(x, y,
         color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=2,
         label='current model')

# ------ format y axis ------ #
yll = y.min()
yhl = y.max()
axes.set_ylim([yll, yhl])
# axes.set_yticks(np.linspace(yll, yhl, 5))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')


# ------ format x axis ------ #
plt.xlim(x.min(), x.max())

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'Sample plot\n '
          r'$x_{range}$' + f'={x_range}'
          f'; \t'
          r'$x_{step}$' + f'={step:.4f}')
plt.xlabel(r'$x_{label}$', fontsize=22)
plt.ylabel(r'$y_{label}$', fontsize=22)
plt.legend()
plt.grid()

# Major ticks every 20, minor ticks every 5
major_ticks = np.arange(0, x_range, 20*step)
minor_ticks = np.arange(0, x_range, 5*step)

axes.set_xticks(major_ticks)
axes.set_xticks(minor_ticks, minor=True)
axes.set_yticks(major_ticks)
axes.set_yticks(minor_ticks, minor=True)

# And a corresponding grid
axes.grid(which='both')

# Or if you want different settings for the grids:
axes.grid(which='minor', alpha=0.2)
axes.grid(which='major', alpha=0.5)
# plt.grid(True)

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure

