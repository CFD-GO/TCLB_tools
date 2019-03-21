import os
import sys
import pandas as pd
import unittest

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir
from DataIO.VTIFile import VTIFile


class TestVtiReader(unittest.TestCase):

    def test_vtk_reader_raise_file_not_found(self):
        filename = 'ghost_file.vti'
        with self.assertRaises(FileNotFoundError) as context:
            vti_reader = VTIFile(filename)

    def test_path(self):
        wd = os.getcwd()
        raise Exception(f"My WD is {wd}")

    def test_vtk_reader(self):
        # real_path = os.path.realpath(__file__)
        # dir_path = os.path.dirname(real_path)
        # print(f'dir_path{dir_path}')

        wd = os.getcwd()
        # print(f"wd{wd}")
        # wd = os.path.dirname(wd)  # go level up

        filename = 'laplace_benchmark_d2q9_VTK_P00_00050010.vti'
        if 'circleci' in wd:
            wd = os.path.join(wd, 'tests')

        filepath = os.path.join(wd, 'sample_data_for_vtk_reader', filename)

        print(f"filepath{filepath}")
        vti_reader = VTIFile(filepath)

        T = vti_reader.get("T")
        U = vti_reader.get("U", vector=True)

        # convert to pandas format
        data = T
        pdT = pd.DataFrame(data=data[1:, 1:],  # values
                           index=data[1:, 0],  # 1st column as index
                           columns=data[0, 1:])  # 1st row as the column names

        pdT2 = pd.DataFrame(data)

        assert T.size > 0
        assert U.size > 0
        assert pdT.size > 0
        assert pdT2.size > 0

    def test_txt_reader(self):
        wd = os.getcwd()
        if 'circleci' in wd:
            wd = os.path.join(wd, 'tests')
        # wd = os.path.dirname(wd)  # go levelup

        filepathT = os.path.join(wd, 'sample_data_for_vtk_reader',
                                 'laplace_benchmark_d2q9_TXT_P00_00050010_T.txt')
        filepathU = os.path.join(wd, 'sample_data_for_vtk_reader',
                                 f'laplace_benchmark_d2q9_TXT_P00_00050010_U.txt')

        dataT = pd.read_csv(filepathT, delimiter=" ")
        dataU = pd.read_csv(filepathU, delimiter=" ")

        assert dataT.size > 0
        assert dataU.size > 0
