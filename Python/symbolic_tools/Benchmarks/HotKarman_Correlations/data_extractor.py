import os
import pandas as pd
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein
import numpy as np
import pwd

import re

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'HotBarman3D')

size = 1
Re = 10
Pr = 10
filename_csv = f'HotKarman3D_template_sizer_{size}_Re{Re}_Pr{Pr}_Log_P00_00000000.csv'
filepath_csv = os.path.join(main_folder, f'keepU_sizer_{size}_Re{Re}_Pr{Pr}', filename_csv)

log = pd.read_csv(filepath_csv, delimiter=",")

avg_outlet = log['HeatFluxX'][-100:].mean()
avg_source = log['HeatSource'][-100:].mean()


folders = os.listdir(main_folder)

logs = []
for folder in folders:
    files = os.listdir(os.path.join(main_folder, folder))
    log = [i for i in files if i.endswith('.csv')]
    logs.append(log)

df = pd.DataFrame(columns=['Re', 'Pr', 'Nu_avgLB_source', 'Nu_conv_avg_outlet', 'Nu_corr', 'Nu_matlab'])  # consider adding nu sdt


def calc_Nu(q_conv, k, D, L):
    T_surf = 1
    T_inf = 0
    Surface = np.pi * D * L
    ht_coeff_experimental = q_conv / (Surface * (T_surf - T_inf))
    Nu = ht_coeff_experimental * D / k
    return Nu


for root, dirs, files in os.walk(main_folder):
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

            new_row = pd.DataFrame({'Re': Re, 'Pr': Pr,
                                    'Nu_avgLB_source': Nu_conv_avg_source,
                                    'Nu_conv_avg_outlet': Nu_conv_avg_outlet,
                                    'Nu_corr': Nu_corr},  index=[0])

            # # simply concatenate both dataframes
            # df = pd.concat([new_row, df]).reset_index(drop=True)
            #
            # # for appending df2 at the end of df1
            # df.append(df2, ignore_index=True)

            # print(f"viscosity={v:0.6f} \t tau={tau:0.3f} \t omega={1 / tau:0.3f}")


print("bye")


# todo
# 1 odzyskac Nu(Re,Pr) na podstawie logow
# 2 narysowac Nu_LB i Nu_correlation w funkcji Re. Osobna linia dla roznych liczb Pr


# df.plot.box()
# https://pandas.pydata.org/pandas-docs/stable/user_guide/visualization.html