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

# df = pd.read_pickle("./pickled_df_Da_1.00e+03_dense2sparse_samples_5.pkl")
df = pd.read_pickle("./pickled_df_Da_1.00e-03_dense2sparse_samples_5.pkl")

plot_dir = 'AC_plots_2D_epsx_epst_dense2sparse'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    

fig_name = os.path.join(plot_dir, 
                        f"Da_{eat_dots_for_texmaker(df['Da'][0])}"
                        f"_Pe_{eat_dots_for_texmaker(df['Pe'][0])}"
                        f"_diffusivity0_{eat_dots_for_texmaker(df['M0'][0])}"
                        f"_lambda_ph0_{eat_dots_for_texmaker(df['lambda0'][0])}"
                        )


plt.rcParams.update({'font.size': 20})
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(17,6))

ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]

dx = np.logspace(np.log10(1./df['L'].min()), np.log10(0.5*1./df['L'].max()), 100)
y = dx**1
y = y / y[0] * df['err_L2'].max()
l1 = ax1.loglog(dx,y, label=r'$\mathcal{O}(\Delta x)$', linewidth=2, color='black', linestyle="--")

y = dx**2
y = y / y[0] * df['err_L2'].max()
l2 = ax1.loglog(dx,y, label=r'$\mathcal{O}(\Delta x^2)$', linewidth=2, color='black', linestyle="-")


dt = np.logspace(np.log10(1./df['n_iterations'].min()), np.log10(0.25*1./df['n_iterations'].max()), 100)
y = dt**1
y = y / y[0] * df['err_L2'].max()
ax2.loglog(dt,y, label=r'$\mathcal{O}(\Delta t)$', linewidth=2, color='black', linestyle="--")

dt = np.logspace(np.log10(1./df[(df['scaling'] == 'acoustic_scaling')]['n_iterations'].min()), 
                 1.5* np.log10(1./df[(df['scaling'] == 'acoustic_scaling')]['n_iterations'].max()), 100)
y = dt**2
y = y / y[0] * df['err_L2'].max()
ax2.loglog(dt,y, label=r'$\mathcal{O}(\Delta t^2)$', linewidth=2, color='black', linestyle="-")



dcpu = np.logspace(np.log10(df['CPU_cost'].min()), 1.05*np.log10(df['CPU_cost'].max()),100)
y = dcpu**(2./4.)
y =  y[0]/y * df['err_L2'].max()
l6 = ax3.loglog(dcpu,y, label=r'$\mathcal{O}(n^{1/2})$', linewidth=2, color='black', linestyle="-.")

dcpu = np.logspace(np.log10(df[(df['scaling'] == 'acoustic_scaling')]['CPU_cost'].min()), 
                  1.05*np.log10(df[(df['scaling'] == 'acoustic_scaling')]['CPU_cost'].max()), 100)
y = dcpu**(2./3.)
y = y[0] /y * df['err_L2'].max()
l7 = ax3.loglog(dcpu,y, label=r'$\mathcal{O}(n^{2/3})$', linewidth=2, color='black', linestyle=":")

scalings = ['acoustic_scaling', 'diffusive_scaling']
markers = ['v', 'o']
for scaling, marker in zip(scalings, markers):
    filtered_df = df[
        # (df['Da'] == kernel_filter) &
        (df['scaling'] == scaling)
        ]
                        
    l3 = ax1.loglog(1./filtered_df['L'], filtered_df['err_L2'], linestyle="", color='black', marker=marker, markersize=10)
    l4 = ax2.loglog(1./filtered_df['n_iterations'], filtered_df['err_L2'],  linestyle="", color='black', marker=marker, markersize=10)
    l5 = ax3.loglog(filtered_df['CPU_cost'], filtered_df['err_L2'],  linestyle="", color='black', marker=marker, markersize=10, label=scaling)


ax1.set_xlim([1E-3, 1E-1])
ax1.set_xticks([1E-3, 1E-2, 1E-1])
ax1.set_ylim([1E-4,1E-1])
# ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.set(xlabel=r'$\Delta x$', ylabel=r'$\mathcal{L}_2 \left(\phi(\Delta x, \Delta t), \phi_{fine} \right)$')
ax1.legend(loc='lower left', shadow=False,  bbox_to_anchor=(0.525,0.01))
ax1.grid(which='major')

ax2.set_xticks([1E-5, 1E-4, 1E-3, 1E-2])
ax2.set_xlim([1E-5, 1E-2])
ax2.set_ylim([1E-4,1E-1])
ax2.set(xlabel=r'$\Delta t$') #, ylabel=r'$L_2(\phi(\delta t), \phi(\delta t_{min})$')
ax2.legend(loc='lower left', shadow=False,  bbox_to_anchor=(0.525,0.01))
ax2.grid(which='major')
ax2.set_yticklabels([]) # Turn off tick labels


ax3.set_xlim([1E5, 1E10])
ax3.set_xticks([1E5, 1E6, 1E7, 1E8, 1E9, 1E10])
ax3.set_ylim([1E-4,1E-1])
ax3.set(xlabel=r'$CPU\; cost$')#, ylabel=r'$L_2(\phi(\delta t), \phi(\delta t_{min})$')
ax3.legend(loc='lower left', shadow=False, bbox_to_anchor=(0.001,0.01))
ax3.grid(which='major')
ax3.set_yticklabels([]) # Turn off tick labels

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.pause(1e-9)  # there is a race condition somewhere in the matplotlib code.
fig.savefig(fig_name + '.pdf', bbox_inches='tight', dpi=200)
fig.savefig(fig_name + '.png', bbox_inches='tight', dpi=200) # for preview
# plt.show()
plt.close(fig)  # close the figure