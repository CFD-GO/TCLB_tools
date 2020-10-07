from DataIO.VTIFile import VTIFile
import os
import pandas as pd

wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

filename = 'laplace_benchmark_d2q9_VTK_P00_00050010.vti'
filepath = os.path.join(wd, 'tests', 'sample_data_for_vtk_reader', filename)
vti_reader = VTIFile(filepath)

T = vti_reader.get("T")
(Ux, Uy, Uz) = vti_reader.get("U", is_vector=True)

# convert to pandas format
data = T
pdT = pd.DataFrame(data=data[1:, 1:],  # values
                   index=data[1:, 0],  # 1st column as index
                   columns=data[0, 1:])  # 1st row as the column names

pdT2 = pd.DataFrame(data)

print('bye')
