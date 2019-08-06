import os
import pandas as pd
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein
import numpy as np
import pwd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import re

#######################################################
#                       README
# to see OS environment variables from python script,
# you need to put them in ~/.profile
#
#######################################################

home = pwd.getpwuid(os.getuid()).pw_dir
local_logs_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'logs_pro')
host_folder = os.path.join(f'$FROM_PRO', 'batch_HotKarman3D_CM_HIGHER', 'keep_nu_and_k_sizer*')

if not os.path.exists(local_logs_folder):
    os.makedirs(local_logs_folder)

cmd = "rsync -zarv  --prune-empty-dirs --include \"*/\"  --include=\"*.csv\" --exclude=\"*\" "\
      + f"\"{host_folder}\""\
      + f" \"{local_logs_folder}\""

print(cmd)
os.system(cmd)


# folders = os.listdir(local_logs_folder)
#
# logs = []
# for folder in folders:
#     files = os.listdir(os.path.join(local_logs_folder, folder))
#     log = [i for i in files if i.endswith('.csv')]
#     logs.append(log)


def calc_Nu(q_conv, k, D, L):
    T_surf = 1
    T_inf = 0
    Surface = np.pi * D * L
    ht_coeff_experimental = q_conv / (Surface * (T_surf - T_inf))
    Nu = ht_coeff_experimental * D / k
    return Nu


def make_plot(x, y, x2, y2, fig_name):
    if not os.path.exists('convergence_plots'):
        os.makedirs('convergence_plots')

    params = {'legend.fontsize': 'xx-large',
              'figure.figsize': (14, 8),
              'axes.labelsize': 'xx-large',
              'axes.titlesize': 'xx-large',
              'xtick.labelsize': 'xx-large',
              'ytick.labelsize': 'xx-large'}
    pylab.rcParams.update(params)

    axes = plt.gca()
    plt.plot(x, y,
             color="black", marker="", markevery=25, markersize=7, linestyle="-", linewidth=2,
             label='Heater')

    plt.plot(x2, y2,
             color="black", marker=">", markevery=25, markersize=7, linestyle=":", linewidth=2,
             label='Outlet')
    # ------ format y axis ------ #
    yll = min(y.min(), y2.min())
    yhl = max(y.max(), y2.max())
    # axes.set_ylim([yll, yhl])

    # axes.set_yticks(np.linspace(yll, yhl, 5))
    # axes.set_yticks(np.arange(yll, yhl, 1E-2))
    # axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
    # axes.yaxis.set_major_formatter(xfmt)

    plt.yscale('log')

    # ------ format x axis ------ #
    # plt.xlim(0, x.max())
    # plt.xlim(0, 4e6)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    plt.title(f'')
    plt.xlabel(r'$iterations$')
    plt.ylabel(r'$|\frac{\sigma}{\mu}|$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure


print("--- rsync complete, lets plot it! ---")

for root, dirs, files in os.walk(local_logs_folder):
    for file in files:
        if file.endswith('.csv') and 'toolbox' not in root:
            # print(file)
            filepath = os.path.join(root, file)
            log = pd.read_csv(filepath, delimiter=",")

            # get Re, Pr numbers from LB log
            match = re.search('sizer_(\d[\.\d]?\d?)_', file, re.IGNORECASE)
            size = float(match.group(1))

            u = 0.01/size
            D = 30*size
            v = log['nu'][0]
            Re = u*D/v

            k = log['conductivity-DefaultZone'][0]
            rho = 1
            cp = 1
            Pr = v*rho*cp/k

            # calculate Nu in all possible ways
            q_conv_avg_source = log['HeatSource'][-100:].mean()
            q_conv_avg_outlet = log['HeatFluxX'][-100:].mean()

            L = 3
            Nu_conv_avg_source = calc_Nu(q_conv_avg_source, k, D, L)
            Nu_conv_avg_outlet = calc_Nu(q_conv_avg_outlet, k, D, L)

            Nu_corr = get_Nu_cylinder_by_Churchill_Bernstein(Re=Re, Pr=Pr)

            print(f"\n\n Re={Re:0.1f} Pr={Pr:0.1f} size={size} "
                  f"Nu_conv_avg_source={Nu_conv_avg_source:0.2f} Nu_conv_avg_outlet={Nu_conv_avg_outlet:0.2f} Nu_corr={Nu_corr:0.2f} "
                  )

            def calc_norm_std(x):
                result = abs(np.std(x) / np.mean(x))
                return result

            window = 50
            skip_first_iterations = 100
            # y1 = log['HeatSource'][100:].rolling(window).apply(lambda x: np.std(x) / np.mean(x)).dropna()
            y1 = log['HeatSource'][skip_first_iterations:].rolling(window).apply(calc_norm_std, raw=False).dropna()
            # x1 = np.linspace(start=window - 1, stop=len(y1) + 1, num=len(y1), endpoint=True)

            y2 = log['HeatFluxX'][skip_first_iterations:].rolling(window).apply(calc_norm_std, raw=True).dropna()
            # x2 = np.linspace(start=window - 1, stop=len(y2) + 1, num=len(y2), endpoint=True)

            x = log['Iteration'][skip_first_iterations + window - 1:]

            plot_name = "convergence_plots/convergence_MA_" + re.sub(r".csv", r".png", file)
            plot_name = re.sub(r"_Log_P00_00000000", r"", plot_name)
            make_plot(x, y1, x, y2, plot_name)

            print(f"processed {file}.")


print("bye")
