import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import os
# -------------------- prepare dummy data --------------------

# square fun version
x1 = 0.5
x2 = 31.5
xm = 0.5*(x1+x2)
ym = 1
a = ym / ((xm - x1)*(xm-x2))
step = 0.01
x = np.arange(x1, x2, step)
y = a*(x - x1)*(x-x2)

# sin fun version
A = np.pi/(x2-x1)
B = x1
y2 = np.sin(A*x - A*B)
# y = -4 * x * (x - x_range) / (x_range * x_range)

if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/sample_plot2D_param_b={x2-x1}.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()
plt.plot(x, y,
         color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=2,
         label='current model')

plt.plot(x, y2,
         color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
         label='another model')
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
plt.xlim(x1-0.5, x2+0.5)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'Sample plot\n '
          r'$x_{1}$=' + f'{x1}' + '\t' + r'$x_{2}$=' + f'{x2}'
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

