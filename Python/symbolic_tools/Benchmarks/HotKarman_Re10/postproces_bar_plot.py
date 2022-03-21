import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.pylab as pylab
import os

# df_Pr10 = pd.read_csv('df_Pr10.csv')
# df_Pr100 = pd.read_csv('df_Pr100.csv')
# df_Pr1000 = pd.read_csv('df_Pr1000.csv')

new_df = pd.read_csv('EnhancedTable.csv')
legend_df = pd.read_csv('EnhancedTableLegend.csv')



df = new_df.transpose()
case_sizes = df.iloc[0]
df.rename(columns=case_sizes, inplace=True)
df.drop(df.index[0], inplace=True)
case_sizes = list(case_sizes)

# case_sizes_set = case_sizes[:3]     # Out[13]: ['Pr10_small', 'Pr10_medium', 'Pr10_large']
# yminmax= [0, max(df[case_sizes_set[0]])]
# fig_name = f'plots/kernel_bar_plot_Pr10'

# case_sizes_set = case_sizes[3:6]    # Out[10]: ['Pr100_small', 'Pr100_medium', 'Pr100_large']
# yminmax= [0, max(df[case_sizes_set[0]])]
# fig_name = f'plots/kernel_bar_plot_Pr100'

case_sizes_set = case_sizes[6:]     # Out[12]: ['Pr1000_small', 'Pr1000_medium', 'Pr1000_large']
yminmax= [min(df[case_sizes_set[0]]), max(df[case_sizes_set[0]])]
fig_name = f'plots/kernel_bar_plot_Pr1000'



case_names = list(new_df.columns.values[1:])
case_names = [w.replace('Nu_CM_HIGHER_1st_order_bc', '$CM-TRT_{1st\; order\; bc}$') for w in case_names]
case_names = [w.replace('Nu_CM_HIGHER_2nd_order_bc', '$CM-TRT_{2nd\; order\; bc}$') for w in case_names]

case_names = [w.replace('Nu_Cumulants_1st_order_bc', '$Cumulants-1^{st}_{1st\; order\; bc}$') for w in case_names]
case_names = [w.replace('Nu_Cumulants_2nd_order_bc', '$Cumulants-1^{st}_{2nd\; order\; bc}$') for w in case_names]

case_names = [w.replace('Nu_CM_1st_order_bc', '$CM-1^{st}_{1st\; order\; bc}$') for w in case_names]
case_names = [w.replace('Nu_CM_2nd_order_bc', '$CM-1^{st}_{2nd\; order\; bc}$') for w in case_names]

case_names = [w.replace('Nu_BGK_1st_order_bc', '$CM-SRT_{1st\; order\; bc}$') for w in case_names]
case_names = [w.replace('Nu_BGK_2nd_order_bc', '$CM-SRT_{2nd\; order\; bc}$') for w in case_names]




# ------------------------------------ PLOT ------------------------------------
if not os.path.exists('plots'):
    os.makedirs('plots')


params = {'legend.fontsize': 'xx-large',
          # 'figure.figsize': (14, 8),
        'figure.figsize': (18, 9),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
#

bar_width = 0.75


pos = np.arange(len(case_names))
plt.subplot(1,3,1)
for i in range(len(case_names)):
    plt.bar(x=pos[i], height=df[case_sizes_set[0]][i], width=bar_width, label=case_names[i])
plt.xticks([],[])
plt.ylim(yminmax)
plt.ylabel('Nu')
plt.grid(True,  axis='y')
plt.title(case_sizes_set[0], fontsize=22)

#The below code will create the second plot.
plt.subplot(1,3,2)
for i in range(len(case_names)):
    plt.bar(x=pos[i], height=df[case_sizes_set[1]][i], width=bar_width, label=case_names[i])
plt.xticks([], [])
plt.ylim(yminmax)
plt.grid(True,  axis='y')
plt.title(case_sizes_set[1], fontsize=22)

#The below code will create the third plot.
ax = plt.subplot(1,3,3)
for i in range(len(case_names)):
    plt.bar(x=pos[i], height=df[case_sizes_set[2]][i], width=bar_width, label=case_names[i])
plt.xticks([],[])
plt.ylim(yminmax)
plt.grid(True,  axis='y')
plt.title(case_sizes_set[2], fontsize=22)


fig = plt.gcf()  # get current figure
# fig.suptitle('Influence of kernel and BC on Nu number', fontsize=26)
fig.tight_layout()
fig.subplots_adjust(bottom=0.2)
# https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
ax.legend(loc='upper center', bbox_to_anchor=(-0.725, -0.05),
          fancybox=True, shadow=True, ncol=5)
plt.show()
# plt.xlabel(r'$edge \, length \, [lu]}$', fontsize=22)
# plt.ylabel(r'$mem \, to \, allocate\, [GB]$', fontsize=22)
#
fig.savefig(fig_name+'.png', bbox_inches='tight')
fig.savefig(fig_name+'.pdf', bbox_inches='tight')
plt.show()
