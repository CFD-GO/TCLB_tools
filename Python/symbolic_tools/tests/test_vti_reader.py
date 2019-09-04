import os
import sys
import pandas as pd
import unittest
import unittest
import matplotlib.pyplot as plt
import numpy as np

import hashlib

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir
from DataIO.VTIFile import VTIFile


class TestVtiReader(unittest.TestCase):

    def test_vtk_reader_raise_file_not_found(self):
        filename = 'ghost_file.vti'
        with self.assertRaises(FileNotFoundError) as context:
            vti_reader = VTIFile(filename)

    def test_vtk_reader_2D(self):
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
        U = vti_reader.get("U", is_vector=True)

        # convert to pandas format
        data = T
        pdT = pd.DataFrame(data=data[1:, 1:],  # values
                           index=data[1:, 0],  # 1st column as index
                           columns=data[0, 1:])  # 1st row as the column names

        pdT2 = pd.DataFrame(data)

        assert T.size > 0
        assert U[0].size > 0
        assert U[1].size > 0
        assert U[2].size > 0
        assert pdT.size > 0
        assert pdT2.size > 0
        filepath = os.path.join(wd, 'sample_data_for_vtk_reader', 'laplace_benchmark_d2q9_TXT_P00_00050010_T.txt')
        T_read_by_pd = pd.read_csv(filepath, delimiter=" ", header=None)
        np.testing.assert_allclose(pdT2 - T_read_by_pd, 0, rtol=1e-14, atol=1e-14)

    def test_vtk_parallel_reader_3D(self):
        filename = 'HotKarman3D_template_sizer_1_nu_0.03_k_0.003_VTK_P00_00600000.pvti'
        wd = os.getcwd()
        # wd = os.path.dirname(wd)  # go level up

        if 'circleci' in wd:
            wd = os.path.join(wd, 'tests')
        filepath = os.path.join(wd, 'sample_data_for_vtk_reader', '3D_multiGPU', filename)

        vti_reader = VTIFile(filepath, parallel=True)
        T_num = vti_reader.get("T")
        [ux_num, uy_num, uz_num] = vti_reader.get("U", is_vector=True)

        expected_hashes = ['5849cd96453a0452c80e37d582fca19f',
                           '20f827ffad4ad10aa50839797316a0eb',
                           '3533e058238c125fcf00592c2269e3d4',
                           '25436c02e5c33da5c8da71338248c423']

        zzs = [T_num, ux_num, uy_num, uz_num]

        hasher = hashlib.md5()
        for expected_hash, zz in zip(expected_hashes, zzs):
            hasher.update(zz.copy(order='C'))
            hash = hasher.hexdigest()
            assert expected_hash == hash


    def test_txt_reader(self):
        wd = os.getcwd()
        if 'circleci' in wd:
            wd = os.path.join(wd, 'tests')
        # wd = os.path.dirname(wd)  # go levelup

        filepathT = os.path.join(wd, 'sample_data_for_vtk_reader',
                                 'laplace_benchmark_d2q9_TXT_P00_00050010_T.txt')
        filepathU = os.path.join(wd, 'sample_data_for_vtk_reader',
                                 f'laplace_benchmark_d2q9_TXT_P00_00050010_U.txt')

        dataT = pd.read_csv(filepathT, delimiter=" ", header=None)
        dataU = pd.read_csv(filepathU, delimiter=" ", header=None)

        assert dataT.size > 0
        assert dataU.size > 0
