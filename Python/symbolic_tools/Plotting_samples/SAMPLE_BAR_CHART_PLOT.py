import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.pylab as pylab
import os

edge_length = np.array([32, 64, 128, 256, 512, 1024])
no_nodes = edge_length*edge_length*edge_length

mem_per_node_B = 976  # Bytes
GB_per_double = 8E-9

mem_to_allocate_GB = no_nodes*mem_per_node_B*GB_per_double/8


description = [str(i) for i in edge_length]
index = np.arange(len(description))

# ------------------------------------ PLOT ------------------------------------
if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/sample_bar_log_plot.pdf'

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (14, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

axes = plt.gca()

plt.bar(index, mem_to_allocate_GB, color="darkblue")
plt.xticks(index, description)
axes.set_yscale('log')

plt.xlabel(r'$edge \, length \, [lu]}$', fontsize=22)
plt.ylabel(r'$mem \, to \, allocate\, [GB]$', fontsize=22)

plt.grid(True,  axis='y')
fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# a = np.array([5, 10, 50, 100, 500, 1000])
#
# p = [0.3, 0.5, 0.2]
# c = np.c_[p[0] * a, p[1] * a, p[2] * a]

# d = np.zeros(c.shape)
# for j, row in enumerate(c):
#     g = np.zeros(len(row) + 1)
#     G = np.sum(row)
#     g[1:] = np.cumsum(row)
#     f = 10 ** (g / G * np.log10(G))
#     f[0] = 0
#     d[j, :] = np.diff(f)
#
# collabels = ["{:3d}%".format(int(100 * i)) for i in p]
#
# dfo = pd.DataFrame(c, columns=collabels)
# df2 = pd.DataFrame(d, columns=collabels)
#
# fig, axes = plt.subplots(ncols=2)
#
# axes[0].set_title("linear stack bar")
# dfo.plot.bar(stacked=True, log=False, ax=axes[0])
# axes[0].set_xticklabels(a)

# axes[1].set_title("log total barheight\nlinear stack distribution")
# df2.plot.bar(stacked=True, log=True, ax=axes[1])
# axes[1].set_xticklabels(a)
# axes[1].set_ylim([1, 1100])
# plt.show()
