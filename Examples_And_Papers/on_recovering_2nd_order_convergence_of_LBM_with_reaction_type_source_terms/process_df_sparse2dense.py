#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 18:01:26 2021

@author: grzegorz
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker
import os, sys
import numpy as np
import re

def eat_dots_for_texmaker(value):
    # s_value = str(value)
    s_value = f"{value:.2e}"
    s_value = re.sub(r"\.", 'o', s_value)
    return s_value

df = pd.read_pickle("./pickled_df_Da_1.00e+03_sparse2dense_samples_4.pkl")

plot_dir = 'AC_plots_2D_epsx_epst'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    

fig_name = os.path.join(plot_dir, 
                        f"Da_{eat_dots_for_texmaker(df['Da'][0])}"
                        f"_Pe_{eat_dots_for_texmaker(df['Pe'][0])}"
                        f"_diffusivity0_{eat_dots_for_texmaker(df['M0'][0])}"
                        f"_lambda_ph0_{eat_dots_for_texmaker(df['lambda0'][0])}"
                        )


plt.rcParams.update({'font.size': 18})


# fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))    
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

# ax1 = axs
ax1 = axs[0]
ax2 = axs[1]

# dx = np.logspace(np.log10(df['LBMdx'].max()), np.log10(df['LBMdx'].min()),100)
# dt = np.logspace(np.log10(df['LBMdt'].max()), np.log10(df['LBMdt'].min()),100)

dx = np.logspace(np.log10(df['L'].min()), np.log10(df['L'].max()), 100)
dt = np.logspace(np.log10(df['n_iterations'].min()), np.log10(df['n_iterations'].max()),100)



y = dx**1
# y = y / y[0] * df['err_L2'].max()
y = y[0]/y * df['err_L2'].max()
l1 = ax1.loglog(dx,y, label=r'${x}$', linewidth=2, color='black', linestyle="--")


y = dx**2
y = y / y[0] * df['err_L2'].max()
y = y[0]/y * df['err_L2'].max()
l2 = ax1.loglog(dx,y, label=r'${x^2}$', linewidth=2, color='black', linestyle="-")

y = dt**1
# y = y / y[0] * df['err_L2'].max()
y =  y[0]/y * df['err_L2'].max()
ax2.loglog(dt,y, label=r'${t}$', linewidth=2, color='black', linestyle="--")

y = dt**2
# y = y / y[0] * df['err_L2'].max()
y = y[0] /y * df['err_L2'].max()
ax2.loglog(dt,y, label=r'${t^2}$', linewidth=2, color='black', linestyle="-")


scalings = ['acoustic_scaling', 'diffusive_scaling']
markers = ['v', 'o']
for scaling, marker in zip(scalings, markers):
    filtered_df = df[
        # (df['Da'] == kernel_filter) &
        (df['scaling'] == scaling)
        ]
                        
    # ax1.loglog(filtered_df['LBMdx'], filtered_df['err_L2'], linestyle="", color='black', marker=marker, markersize=10, label=scaling)
    # ax2.loglog(filtered_df['LBMdt'], filtered_df['err_L2'],  linestyle="", color='black', marker=marker, markersize=10, label=scaling)
    
    l3 = ax1.loglog(filtered_df['L'], filtered_df['err_L2'], linestyle="", color='black', marker=marker, markersize=10)
    l4 = ax2.loglog(filtered_df['n_iterations'], filtered_df['err_L2'],  linestyle="", color='black', marker=marker, markersize=10)


ax1.set_xticks([1E1, 1E2, 1E3])
ax2.set_ylim([1E-4,2E-1])
ax1.set_yticks([1E-1, 1E-2, 1E-3, 1E-4])
# ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.set(xlabel=r'$L$', ylabel=r'$L_2(\phi(\delta x), \phi(\delta x_{min})$')
# ax1.set(xlabel=r'$\epsilon_x$', ylabel=r'$L_2(\phi(\delta x), \phi(\delta x_{min})$')
# ax1.set_xscale('log', base=2)
# ax1.legend()
ax1.grid(which='both')
# ax1.grid(True)

ax2.set_xticks([5E2, 1E4, 5E4])
ax2.set_xticks([1E3, 1E4])
ax2.set_ylim([1E-4,1E-1])
ax2.set_yticks([1E-1, 1E-2, 1E-3, 1E-4])
ax2.set(xlabel=r'$iterations$', ylabel=r'$L_2(\phi(\delta t), \phi(\delta t_{min})$')
# ax2.set(xlabel=r'$\epsilon_t$', ylabel=r'$L_2(\phi(\delta t), \phi(\delta t_{min})$')
# ax2.set_xscale('log', base=2)
# ax2.legend()
ax2.grid(which='both')
# ax2.grid(True)

line_labels = ['$\mathcal{O}(n)$', '$\mathcal{O}(n^2)$', "acoustic scaling", "diffusive scaling"]
# Create the legend
fig.legend([l1, l2, l3, l4],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="lower left",   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(0.75,0.65)
           # title="Legend Title"  # Title for the legend
           )

# Adjust the scaling factor to fit your legend text completely outside the plot
# (smaller value results in more space being made for the legend)
# plt.subplots_adjust(right=0.8)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.pause(1e-9)  # there is a race condition somewhere in the matplotlib code.
fig.savefig(fig_name + '.pdf', bbox_inches='tight', dpi=200)
fig.savefig(fig_name + '.png', bbox_inches='tight', dpi=200) # for preview
# plt.show()
plt.close(fig)  # close the figure