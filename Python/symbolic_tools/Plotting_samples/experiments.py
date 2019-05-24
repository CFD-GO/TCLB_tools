import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# -------------------- prepare dummy data --------------------

x_range = 1
step = 0.01
x = np.arange(0.0, x_range, step)
# y = -4 * x * (x - x_range) / (x_range * x_range)
fig_name = f'sample_plot2D_param_b={x_range}.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()
plt.plot(x, 1/(x+1),
         color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=2,
         label=r'$f^{wall}_{\bar{i}}$')

plt.plot(x, x/(x+1),
         color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
         label=r'$f^{\star}_{\bar{i}}$')
# plt.plot(x, 1/(2*x),
#          color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=2,
#          label=r'$f^{\star}_i$')
#
# plt.plot(x, (2*x-1)/(2*x),
#          color="black", marker="", markevery=1, markersize=15, linestyle="-", linewidth=2,
#          label=r'$f^{\star}_{\bar{i}}$')

# ------ format y axis ------ #
# yll = y.min()
# yhl = y.max()
# axes.set_ylim([yll, yhl])
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
plt.xlabel(r'$x_{label}$')
plt.ylabel(r'$y_{label}$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure

