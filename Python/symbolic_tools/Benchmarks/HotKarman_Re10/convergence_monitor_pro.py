import os
import pandas as pd
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein
import numpy as np
import pwd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import re
import warnings
from fractions import Fraction


#######################################################
#                       README
# to see OS environment variables from python script,
# you need to put them in ~/.profile
#
#######################################################

home = pwd.getpwuid(os.getuid()).pw_dir
local_logs_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotKarman_Re10')
# host_folder = os.path.join(f'$FROM_PRO', 'batch_HotKarman_Re10_test', 'batch_HotKarman_Re10_sizer*')
#
# if not os.path.exists(local_logs_folder):
#     os.makedirs(local_logs_folder)
#
# cmd = "rsync -zarv  --prune-empty-dirs --include \"*/\"  --include=\"*.csv\" --exclude=\"*\" "\
#       + f"\"{host_folder}\""\
#       + f" \"{local_logs_folder}\""
#
# print(cmd)
# os.system(cmd)
print("--- rsync complete ---")


# folders = os.listdir(local_logs_folder)
#
# logs = []
# for folder in folders:
#     files = os.listdir(os.path.join(local_logs_folder, folder))
#     log = [i for i in files if i.endswith('.csv')]
#     logs.append(log)


def calc_Nu(q_conv, k, D, L):
    T_surf = 11
    T_inf = 10
    Surface = np.pi * D * L
    ht_coeff_experimental = q_conv / (Surface * (T_surf - T_inf))
    Nu = ht_coeff_experimental * D / k
    # Nu_jfm = q_conv*D/(k*L*(T_surf - T_inf))
    return Nu


def make_std_q_plot(x, y, x2, y2, fig_name):
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


def make_Nu_plot(x, y, x2, y2, fig_name, plot_title=''):
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
             label='Source')

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

    # plt.yscale('log')

    # ------ format x axis ------ #
    # plt.xlim(0, x.max())
    # plt.xlim(0, 4e6)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)

    plt.title(f'{plot_title}')
    plt.xlabel(r'$iterations$')
    plt.ylabel(r'$Nu$')
    plt.legend()
    plt.grid()

    fig = plt.gcf()  # get current figure
    fig.savefig(fig_name, bbox_inches='tight')
    plt.show()

    plt.close(fig)  # close the figure


def get_data_from_path(filepath):
    log = pd.read_csv(filepath, delimiter=",")

    # get Re, Pr numbers from LB log
    match = re.search(r'Re(\d{1,4})_', file, re.IGNORECASE)
    Re = float(match.group(1))

    match = re.search(r'_D0_(\d\.\d\de\+\d\d)_', file, re.IGNORECASE)
    D0 = float(match.group(1))

    match = re.search(r'_U_(\d\.\d\de\-\d\d)_', file, re.IGNORECASE)
    u = float(match.group(1))

    match = re.search(r'Pr_?(\d{1,4})_', file, re.IGNORECASE)
    Pr = float(match.group(1))

    match = re.search(r'_BR_(\d\.\d\de\-\d\d)_', file, re.IGNORECASE)
    if match is not None:
        blockage_ratio = float(match.group(1))
        blockage_ratio = Fraction(blockage_ratio).limit_denominator(100)
        blockage_ratio = str(blockage_ratio).replace('/', 'over')
    match = re.search(r'_BR_(\do\d{1,4})_', file, re.IGNORECASE)
    if match is not None:
        blockage_ratio = match.group(1)
    else:
        blockage_ratio = 'unknownBR'  # d_cylinder / domain_y

    if bool(re.search("is3d", filepath)):
        is3d = bool((re.search("is3d_TRUE", filepath)))
    else:
        raise Exception("is it 2d or 3d simulation?")

    match = re.search(r'\d([a-z]){2}_order_bc', file, re.IGNORECASE)
    bc_order = match.group(0)

    match = re.search(r'Cumulants|CM_HIGHER', file, re.IGNORECASE)
    collision_kernel = match.group(0)

    re.search(r'Cumulants', file, re.IGNORECASE)

    v = log['nu'][0]
    # if Re != u*D0/v:
    #     print(f"Re={Re:0.1f} while u*D0/v = {u*D0/v:0.1f}")

    k = log['conductivity-DefaultZone'][0]
    rho = 1
    cp = 1
    np.testing.assert_allclose(Pr - v * rho * cp / k, 0, rtol=1e-2, atol=1e-2)

    # check for nan
    nanstr = ' -nan'
    # mask = log['HeatFluxX2'].values == ' -nan'
    mask = (log['HeatSource'].values != nanstr) & (log['HeatFluxX'].values != nanstr) & (log['HeatFluxX2'].values != nanstr)

    if np.any(mask) and isinstance(mask, np.ndarray):
        warnings.warn("NaN detected. For further analysis NaN has been sliced out from log", UserWarning)
        log = log[mask].astype('float')

    df = pd.DataFrame({'Re': Re, 'Pr': Pr, 'BR': blockage_ratio, 'D': D0, 'nu': v, 'U': u, 'Collision_Kernel': collision_kernel, 'BC_order': bc_order, 'is3D': is3d, 'log_length': len(log)}, index=[0])
    return log, df


print("--- lets plot it! ---")
# dfObj = pd.DataFrame(columns=['Re', 'Pr', 'BR', 'nu', 'U', 'BC_order', 'Collision_Kernel', 'is3D'])
Nu_from_QS_FEM = {
    'Pr10': 4.82,
    'Pr100': 10.10,
    'Pr1000': 21.43,
}

dfObj = pd.DataFrame()

for root, dirs, files in os.walk(local_logs_folder):
    for file in files:
        if file.endswith('.csv') and 'toolbox' not in root:
            # print(file)
            filepath = os.path.join(root, file)
            log, df = get_data_from_path(filepath)


            # calculate Nu in all possible ways
            n_last_it_to_analyze = 1000
            q_conv_avg_source = log['HeatSource'][-n_last_it_to_analyze:].mean()
            q_conv_avg_inlet = log['HeatFluxX'][-n_last_it_to_analyze:].mean()
            q_conv_avg_outlet = log['HeatFluxX2'][-n_last_it_to_analyze:].mean()
            dq_conv_avg_outlet = -(q_conv_avg_outlet - q_conv_avg_inlet)

            L = (1, 3)[df['is3D'][0]]
            k = df['nu'][0]/df['Pr'][0]  # assuming rho=cp=1
            Nu_conv_avg_source = calc_Nu(q_conv_avg_source, k, df['D'][0], L)
            # Nu_conv_avg_inlet = calc_Nu(q_conv_avg_inlet, k, D0, L)
            # Nu_conv_avg_outlet = calc_Nu(q_conv_avg_outlet, k, D0, L)
            Nu_conv_avg_doutlet = calc_Nu(dq_conv_avg_outlet, k, df['D'][0], L)

            Nu_corr = get_Nu_cylinder_by_Churchill_Bernstein(Re=df['Re'][0], Pr=df['Pr'][0])
            Pr_str = "".join(('Pr', str(int(df['Pr']))))
            Nu_txt_description = f"Nu_avg_source={Nu_conv_avg_source:0.2f} " \
                                 f"Nu_avg_outlet={Nu_conv_avg_doutlet:0.2f} " \
                                 f"Nu_FEM={Nu_from_QS_FEM[Pr_str]}"

            print(f"\n\nRe={df['Re'][0]:0.1f} Pr={df['Pr'][0]:0.1f} D0={df['D'][0]} \n{Nu_txt_description}")

            y1 = log['HeatSource'][-n_last_it_to_analyze:].apply(calc_Nu, args=(k, df['D'][0], L)).dropna()

            dHeatFluxX = -(log['HeatFluxX2'] - log['HeatFluxX'])
            y2 = dHeatFluxX[-n_last_it_to_analyze:].apply(calc_Nu, args=(k, df['D'][0], L)).dropna()
            x = log['Iteration'][-n_last_it_to_analyze:]
            # plot_name = f"convergence_plots/BR{blockage_ratio}_Nu_" + re.sub(r".csv", r".png", file)
            plot_name = f"convergence_plots/is3d_{df['is3D'][0]}_" + re.sub(r".csv", r".png", file)
            plot_name = re.sub(r"_Log_P\d{2}_\d{8}", r"", plot_name)
            make_Nu_plot(x, y1, x, y2, plot_name, plot_title=Nu_txt_description)

            def calc_norm_std(x):
                result = abs(np.std(x) / np.mean(x))
                return result

            window = 50
            # y1 = log['HeatSource'][100:].rolling(window).apply(lambda x: np.std(x) / np.mean(x)).dropna()
            y1 = log['HeatSource'][-n_last_it_to_analyze:].rolling(window).apply(calc_norm_std, raw=False).dropna()
            # x1 = np.linspace(start=window - 1, stop=len(y1) + 1, num=len(y1), endpoint=True)

            y2 = log['HeatFluxX'][-n_last_it_to_analyze:].rolling(window).apply(calc_norm_std, raw=True).dropna()
            # x2 = np.linspace(start=window - 1, stop=len(y2) + 1, num=len(y2), endpoint=True)

            x = log['Iteration'][-n_last_it_to_analyze + window - 1:]

            plot_name = "convergence_plots/convergence_MA_" + re.sub(r".csv", r".png", file)
            plot_name = re.sub(r"_Log_P{\d\d}_\d{8}", r"", plot_name)
            # make_std_q_plot(x, y1, x, y2, plot_name)

            df = df.assign(Nu_avg_source=round(Nu_conv_avg_source, 2),
                           Nu_avg_outlet=round(Nu_conv_avg_doutlet, 2),
                           Nu_FEM=Nu_from_QS_FEM[Pr_str]
                           )
            dfObj = dfObj.append(df)
            print(f"processed {file}.")


print(dfObj)

dfObj = dfObj.sort_values(by=['is3D', 'Collision_Kernel', 'D', 'BC_order'])
dfObj.to_pickle("./pickled_df.pkl")

# with pd.ExcelWriter('LBM_validation_HotKarman_Benchmark.xlsx') as writer:  # doctest: +SKIP
#     # dfObj = dfObj.sort_values(by=['is3D', 'D', 'BC_order', 'Collision_Kernel', 'log_length'])
#
#     dfObjPr10 = dfObj.loc[dfObj['Pr'] == 10]
#     dfObjPr10.to_excel(writer, sheet_name='Pr10')
#     dfObjPr100 = dfObj.loc[dfObj['Pr'] == 100]
#     dfObjPr100.to_excel(writer, sheet_name='Pr100')

print("bye")
