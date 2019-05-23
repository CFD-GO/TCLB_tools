import os
import pandas as pd
from Benchmarks.HT_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein
import numpy as np

wd = os.getcwd()
wd = os.path.dirname(wd)  # go levelup
wd = os.path.dirname(wd)  # go levelup

# HotKarman2D_template_nu_0.03_k_3e-05_Log_P00_00000000

filepath = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', f'sample_tclb_log.csv')
# log = pd.read_csv(filepath, delimiter=",", header=None)
log = pd.read_csv(filepath, delimiter=",")


avg_outlet = log['HeatFluxX'][-100:].mean()
avg_source = log['HeatSource'][-100:].mean()


path = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', 'batch_HotKarman2D_EQ')
files = os.listdir(path)

logs = [i for i in files if i.endswith('.csv')]
df = pd.DataFrame(columns=['Re', 'Pr', 'Nu_avgLB_source', 'Nu_conv_avg_outlet', 'Nu_corr', 'Nu_matlab'])  # consider adding nu sdt

def calc_Nu(q_conv, k, D, L):
    T_surf = 1
    T_inf = 0
    Surface = np.pi * D * L
    ht_coeff_experimental = q_conv / (Surface * (T_surf - T_inf))
    Nu = ht_coeff_experimental * D / k
    return Nu

for root, dirs, files in os.walk(path):
    for file in files:
        if file.endswith('.csv'):
            print(file)

            log = pd.read_csv(filepath, delimiter=",")
            # get Re, Pr numbers from LB log
            u = 0.01
            D = 30
            v = log['nu'][0]
            Re = u*D/v

            k = log['conductivity-DefaultZone'][0]
            rho = 1
            cp = 1
            Pr = v*rho*cp/k

            # calculate Nu in all possible ways
            q_conv_avg_source = log['HeatSource'][-100:].mean()
            q_conv_avg_outlet = log['HeatFluxX'][-100:].mean()

            L = 3 # get it from case name?
            Nu_conv_avg_source = calc_Nu(q_conv_avg_outlet, k, D, L)
            Nu_conv_avg_outlet = calc_Nu(q_conv_avg_outlet, k, D, L)

            Nu_corr = get_Nu_cylinder_by_Churchill_Bernstein(Re=Re, Pr=Pr)

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