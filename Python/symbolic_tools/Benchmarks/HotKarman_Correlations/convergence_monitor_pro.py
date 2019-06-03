import os
import pandas as pd
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein
import numpy as np
import pwd
import matplotlib.pyplot as plt

import re
home = pwd.getpwuid(os.getuid()).pw_dir
local_logs_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'logs_pro')

USER = "plgmuaddieb"
HOST = "prometheus.cyfronet.pl"

# os.system(f"ssh {USER}@{HOST} ls $SCRATCH/")

print("-------------")
#os.system(f"ssh {USER}@{HOST} ls /net/scratch/people/{USER}")
# main_folder = os.path.join("/net/scratch/people/{USER}", 'output', 'batch_HotKarman3D')


if not os.path.exists(local_logs_folder):
    os.makedirs(local_logs_folder)

cmd = "rsync -zarv  --prune-empty-dirs --include \"*/\"  --include=\"*.csv\" --exclude=\"*\" \"$FROM_PRO/batch_HotKarman3D\"" + f" \"{local_logs_folder}\""

# cmd = f"rsync -zarv  --prune-empty-dirs " \
#     f"--include \"*/\"  " \
#     f"--include=\"*.csv\" " \
#     f"--exclude=\"*\" " \
#     f"\"plgmuaddieb@prometheus.cyfronet.pl:/net/scratch/people/plgmuaddieb/output/batch_HotKarman3D/keep_nuk_old/\" \"{local_logs_folder}\""

print(cmd)
os.system(cmd)


folders = os.listdir(local_logs_folder)
#
logs = []
for folder in folders:
    files = os.listdir(os.path.join(local_logs_folder, folder))
    log = [i for i in files if i.endswith('.csv')]
    logs.append(log)


def calc_Nu(q_conv, k, D, L):
    T_surf = 1
    T_inf = 0
    Surface = np.pi * D * L
    ht_coeff_experimental = q_conv / (Surface * (T_surf - T_inf))
    Nu = ht_coeff_experimental * D / k
    return Nu


def make_plot(x, y, x2, y2, fig_name):
    # -------------------- make dummy plot --------------------
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(14, 8))

    axes = plt.gca()
    plt.plot(x, y,
             color="black", marker="", markevery=12, markersize=7, linestyle="-", linewidth=2,
             label='HeatSource')

    plt.plot(x2, y2,
             color="black", marker=">", markevery=12, markersize=7, linestyle=":", linewidth=2,
             label='HeatFluxX')
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
    plt.xlim(x.min(), x.max())

    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    plt.title(f'Convergence Monitor - Moving Average')
    plt.xlabel(r'$iterations \times 10^3$')
    plt.ylabel(r'$|\frac{\sigma}{\mu}|$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure


for root, dirs, files in os.walk(local_logs_folder):
    for file in files:
        if file.endswith('.csv') and 'toolbox' not in root:
            # print(file)
            filepath = os.path.join(root, file)
            log = pd.read_csv(filepath, delimiter=",")
            # get Re, Pr numbers from LB log
            u = 0.01
            match = re.search('sizer_(\d+)', file, re.IGNORECASE)
            size = int(match.group(1))
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
            Nu_conv_avg_source = calc_Nu(q_conv_avg_outlet, k, D, L)
            Nu_conv_avg_outlet = calc_Nu(q_conv_avg_outlet, k, D, L)

            Nu_corr = get_Nu_cylinder_by_Churchill_Bernstein(Re=Re, Pr=Pr)
            print(f"Re={Re:0.1f} Pr={Pr:0.1f} size={size} "
                  f"Nu_conv_avg_source={Nu_conv_avg_source:0.1f} Nu_conv_avg_outlet={Nu_conv_avg_outlet:0.1f} Nu_corr={Nu_corr:0.1f} "
                  f"--- Extracted from file {file}")

            def calc_norm_std(x):
                result = abs(np.std(x) / np.mean(x))
                return result

            window = 50
            skip_first_iterations = 100
            # y1 = log['HeatSource'][100:].rolling(window).apply(lambda x: np.std(x) / np.mean(x)).dropna()
            y1 = log['HeatSource'][skip_first_iterations:].rolling(window).apply(calc_norm_std).dropna()
            x1 = np.linspace(start=window - 1, stop=len(y1) + 1, num=len(y1), endpoint=True)

            y2 = log['HeatFluxX'][skip_first_iterations:].rolling(window).apply(calc_norm_std).dropna()
            x2 = np.linspace(start=window - 1, stop=len(y2) + 1, num=len(y2), endpoint=True)

            plot_name = "convergence_MA_" + re.sub(r".csv", r".png", file)
            make_plot(x1, y1, x2, y2, plot_name)

            print(f"processed {file}.")
            # y2 = y[window - 1:]
            # x = np.linspace(start=window - 1, stop=len(y2) + 1, num=len(y2), endpoint=True)


print("bye")
