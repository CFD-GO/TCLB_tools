import os
import pandas as pd

wd = os.getcwd()
wd = os.path.dirname(wd)  # go levelup

filepathT = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', 'laplace_benchmark_d2q9_TXT_P00_00050010_T.txt')
filepathU = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', f'laplace_benchmark_d2q9_TXT_P00_00050010_U.txt')

dataT = pd.read_csv(filepathT, delimiter=" ", header=None)
dataU = pd.read_csv(filepathU, delimiter=" ", header=None)

print("bye")
