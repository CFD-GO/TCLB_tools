import unittest

import io
from contextlib import redirect_stdout
from sympy import Symbol

import multiprocessing
from concurrencytest import ConcurrentTestSuite, fork_for_tests

import sys
import os

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, \
    get_mom_vector_from_shift_mat, \
    get_mom_vector_from_discrete_def

from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, \
    F2D, dzeta2D, u2D, \
    rho, w_D2Q9, m00, e_D2Q9, \
    moments_dict

from SymbolicCollisions.core.printers import print_as_vector

from SymbolicCollisions.core.hardcoded_results import \
    hardcoded_F_cm_Guo_hydro_incompressible_D2Q9, hardcoded_cm_eq_compressible_D2Q9


class TestDiscreteCMTransforms(unittest.TestCase):
    def test_shift_vs_def_cm(self):
        dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

        functions = [lambda i: w_D2Q9[i], dcmt.get_force_He, dcmt.get_force_Guo]
        from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, NrawD2Q9

        for fun in functions:
            F_in_cm = get_mom_vector_from_discrete_def(fun,
                                                       discrete_transform=dcmt.get_cm,
                                                       moments_order=moments_dict['D2Q9'],
                                                       serial_run=True)  # calculate from definition of cm
            NMF_cm = get_mom_vector_from_shift_mat(fun, mat=NrawD2Q9 * Mraw_D2Q9)  # calculate using shift matrices

            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(F_in_cm, 'F_in_cm')
            out = f.getvalue()

            f2 = io.StringIO()
            with redirect_stdout(f2):
                print_as_vector(NMF_cm, 'F_in_cm')
            out2 = f2.getvalue()

            assert out == out2

    def test_get_F_cm_Guo_continuous_and_discrete(self):
        dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
        F_cm_Guo_disc = get_mom_vector_from_discrete_def(dcmt.get_force_Guo,
                                                         discrete_transform=dcmt.get_cm,
                                                         moments_order=moments_dict['D2Q9'],
                                                         serial_run=True)

        from SymbolicCollisions.core.ContinuousCMTransforms import \
            ContinuousCMTransforms, get_mom_vector_from_continuous_def

        from SymbolicCollisions.core.cm_symbols import \
            F3D, dzeta3D, u3D

        ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
        F_cm_Guo_cont = get_mom_vector_from_continuous_def(ccmt.get_force_Guo,
                                                           continuous_transformation=ccmt.get_cm,
                                                           moments_order=moments_dict['D2Q9'],
                                                           serial_run=True)

        # print_as_vector(F_cm_Guo_cont, 'F_cm')
        results = [F_cm_Guo_disc, F_cm_Guo_cont]

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(hardcoded_F_cm_Guo_hydro_incompressible_D2Q9, 'F_cm')
        expected_result = f.getvalue()

        for result in results:
            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(result, 'F_cm')
            out = f.getvalue()

            assert out == expected_result

    def test_get_force_He_discrete(self):
        dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
        F_in_cm = get_mom_vector_from_discrete_def(dcmt.get_force_He,
                                                   discrete_transform=dcmt.get_cm,
                                                   moments_order=moments_dict['D2Q9'],
                                                   serial_run=True)
        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_in_cm, 'F_in_cm')
        out = f.getvalue()

        expected_result = '\tF_in_cm[0] = 0;\n' \
                          '\tF_in_cm[1] = F.x*m00/rho;\n' \
                          '\tF_in_cm[2] = F.y*m00/rho;\n' \
                          '\tF_in_cm[3] = -3.*m00*ux2*(F.x*u.x + F.y*u.y)/rho;\n' \
                          '\tF_in_cm[4] = -3.*m00*uy2*(F.x*u.x + F.y*u.y)/rho;\n' \
                          '\tF_in_cm[5] = -3.*m00*uxuy*(F.x*u.x + F.y*u.y)/rho;\n' \
                          '\tF_in_cm[6] = m00*(9.*F.x*ux3*u.y + 9.*F.y*ux2*uy2 + 1/3.*F.y)/rho;\n' \
                          '\tF_in_cm[7] = m00*(9.*F.x*ux2*uy2 + 1/3.*F.x + 9.*F.y*u.x*uy3)/rho;\n' \
                          '\tF_in_cm[8] = -m00*(18.*F.x*ux3*uy2 + F.x*ux3 + 3.*F.x*u.x*uy2 + 18.*F.y*ux2*uy3 + 3.*F.y*ux2*u.y + F.y*uy3)/rho;\n'  # noqa

        assert 'F_in_cm[0] = 0;' in out
        assert 'F_in_cm[1] = F.x*m00/rho;' in out
        assert 'F_in_cm[2] = F.y*m00/rho;' in out
        assert 'F_in_cm[3] = -3.*m00*ux2*(F.x*u.x + F.y*u.y)/rho;\n' in out
        assert 'F_in_cm[4] = -3.*m00*uy2*(F.x*u.x + F.y*u.y)/rho;\n' in out
        assert 'F_in_cm[5] = -3.*m00*uxuy*(F.x*u.x + F.y*u.y)/rho;\n' in out
        assert 'F_in_cm[6] = m00*(9.*F.x*ux3*u.y + 9.*F.y*ux2*uy2 + 1/3.*F.y)/rho;\n' in out
        assert 'F_in_cm[7] = m00*(9.*F.x*ux2*uy2 + 1/3.*F.x + 9.*F.y*u.x*uy3)/rho;\n' in out
        assert 'F_in_cm[8] = -m00*(18.*F.x*ux3*uy2 + F.x*ux3 + 3.*F.x*u.x*uy2 + 18.*F.y*ux2*uy3 + 3.*F.y*ux2*u.y + F.y*uy3)/rho;\n' in out  # noqa

        assert expected_result == out

    def test_get_cm_eq_hydro_discrete(self):
        dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
        cm_eq = get_mom_vector_from_discrete_def(dcmt.get_EDF_incompressible,
                                                 discrete_transform=dcmt.get_cm,
                                                 moments_order=moments_dict['D2Q9'],
                                                 serial_run=True)
        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq')
        out = f.getvalue()

        expected_result = '\tcm_eq[0] = m00;\n' \
                          '\tcm_eq[1] = u.x*(1 - m00);\n' \
                          '\tcm_eq[2] = u.y*(1 - m00);\n' \
                          '\tcm_eq[3] = m00*ux2 + 1/3.*m00 - ux2;\n' \
                          '\tcm_eq[4] = m00*uy2 + 1/3.*m00 - uy2;\n' \
                          '\tcm_eq[5] = uxuy*(m00 - 1.);\n' \
                          '\tcm_eq[6] = u.y*(-m00*ux2 - 1/3.*m00 + 1/3.);\n' \
                          '\tcm_eq[7] = u.x*(-m00*uy2 - 1/3.*m00 + 1/3.);\n' \
                          '\tcm_eq[8] = m00*ux2*uy2 + 1/3.*m00*ux2 + 1/3.*m00*uy2 + 1/9.*m00 + 2.*ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n'  # noqa

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[1] = u.x*(1 - m00);' in out
        assert 'cm_eq[2] = u.y*(1 - m00);' in out
        assert 'cm_eq[3] = m00*ux2 + 1/3.*m00 - ux2;\n' in out
        assert 'cm_eq[4] = m00*uy2 + 1/3.*m00 - uy2;\n' in out
        assert 'cm_eq[5] = uxuy*(m00 - 1.);\n' in out
        assert 'cm_eq[6] = u.y*(-m00*ux2 - 1/3.*m00 + 1/3.);\n' in out
        assert 'cm_eq[7] = u.x*(-m00*uy2 - 1/3.*m00 + 1/3.);\n' in out
        assert 'cm_eq[8] = m00*ux2*uy2 + 1/3.*m00*ux2 + 1/3.*m00*uy2 + 1/9.*m00 + 2.*ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n' in out  # noqa

        assert expected_result == out

    def test_cm_eq_compressible_discrete(self):
        """
        test eq 10 from
        'Modeling incompressible thermal flows using a central-moment-based lattice Boltzmann method'
        Linlin Fei, Kai Hong Luo, Chuandong Lin, Qing Li
        2017
        """
        dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
        cm_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma(i),
                                                 discrete_transform=dcmt.get_cm,
                                                 moments_order=moments_dict['D2Q9'],
                                                 serial_run=True)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq')
        out = f.getvalue()

        expected_result = '\tcm_eq[0] = m00;\n' \
                          '\tcm_eq[1] = 0;\n' \
                          '\tcm_eq[2] = 0;\n' \
                          '\tcm_eq[3] = 1/3.*m00;\n' \
                          '\tcm_eq[4] = 1/3.*m00;\n' \
                          '\tcm_eq[5] = 0;\n' \
                          '\tcm_eq[6] = -m00*ux2*u.y;\n' \
                          '\tcm_eq[7] = -m00*u.x*uy2;\n' \
                          '\tcm_eq[8] = m00*(3.*ux2*uy2 + 1/9.);\n'

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[2] = 0;' in out
        assert 'cm_eq[2] = 0;' in out
        assert 'cm_eq[3] = 1/3.*m00;\n' in out
        assert 'cm_eq[4] = 1/3.*m00;\n' in out
        assert 'cm_eq[5] = 0;\n' in out
        assert 'cm_eq[6] = -m00*ux2*u.y;\n' in out
        assert 'cm_eq[7] = -m00*u.x*uy2;\n' in out
        assert 'cm_eq[8] = m00*(3.*ux2*uy2 + 1/9.);\n' in out

        assert expected_result == out


# Pycharm runs them sequentially
# python -m unittest tests/test_example_unit_tests_parallel_run.py # sequential as well
# python tests/test_example_unit_tests_parallel_run.py # concurrent run :)

if __name__ == '__main__':
    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner()
    # Run same tests across 4 processes
    cores = multiprocessing.cpu_count()
    print(f'\nRunning tests on {cores} cores:')
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDiscreteCMTransforms)
    concurrent_suite = ConcurrentTestSuite(suite, fork_for_tests(cores))
    runner.run(concurrent_suite)
