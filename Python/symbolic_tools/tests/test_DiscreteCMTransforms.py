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
# from Python.symbolic_tools.SymbolicCollisions.core.cm_symbols import w  # alternatively change all import paths

from SymbolicCollisions.core.cm_symbols import w

from SymbolicCollisions.core.DiscreteCMTransforms import \
    get_mom_vector_from_discrete_def, get_mom_vector_from_shift_Mat, \
    get_discrete_EDF_hydro, \
    get_discrete_force_He, \
    get_discrete_force_Guo,\
    get_gamma, get_discrete_cm

from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, NrawD2Q9

from SymbolicCollisions.core.printers import print_as_vector

from SymbolicCollisions.core.hardcoded_results import \
    hardcoded_F_cm_Guo_hydro_LB_velocity_based_D2Q9, hardcoded_cm_eq_compressible_D2Q9


class TestSymbolicCalc(unittest.TestCase):
    def test_shift_vs_def_cm(self):
        functions = [lambda i: w[i], get_discrete_force_He, get_discrete_force_Guo]

        for fun in functions:
            F_in_cm = get_mom_vector_from_discrete_def(fun, discrete_transform=get_discrete_cm)  # calculate from definition of cm
            NMF_cm = get_mom_vector_from_shift_Mat(fun, Mat=NrawD2Q9 * Mraw_D2Q9)  # calculate using shift matrices

            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(F_in_cm, 'F_in_cm', regex=True)
            out = f.getvalue()

            f2 = io.StringIO()
            with redirect_stdout(f2):
                print_as_vector(NMF_cm, 'F_in_cm', regex=True)
            out2 = f2.getvalue()

            assert out == out2


    def test_get_F_cm_Guo_continuous_and_discrete(self):
        F_cm_Guo_disc = get_mom_vector_from_discrete_def(get_discrete_force_Guo, discrete_transform=get_discrete_cm)

        from SymbolicCollisions.core.sym_col_fun import get_mom_vector_from_continuous_def, get_continuous_force_Guo, get_continuous_cm
        F_cm_Guo_cont = get_mom_vector_from_continuous_def(get_continuous_force_Guo, continuous_transformation=get_continuous_cm)

        print_as_vector(F_cm_Guo_cont, 'F_cm', regex=True)

        results = [F_cm_Guo_disc, F_cm_Guo_cont]

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(hardcoded_F_cm_Guo_hydro_LB_velocity_based_D2Q9, 'F_cm', regex=True)
        expected_result = f.getvalue()

        for result in results:
            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(result, 'F_cm', regex=True)
            out = f.getvalue()

        assert out == expected_result

    def test_get_force_He_discrete(self):
        F_in_cm = get_mom_vector_from_discrete_def(get_discrete_force_He, discrete_transform=get_discrete_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_in_cm, 'F_in_cm', regex=True)
        out = f.getvalue()

        expected_result = 'F_in_cm[0] = 0;\n' \
                          'F_in_cm[1] = Fhydro.x*m00/rho;\n' \
                          'F_in_cm[2] = Fhydro.y*m00/rho;\n' \
                          'F_in_cm[3] = -3.*m00*ux2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' \
                          'F_in_cm[4] = -3.*m00*uy2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' \
                          'F_in_cm[5] = -3.*m00*uxuy*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' \
                          'F_in_cm[6] = m00*(9.*Fhydro.x*ux3*u.y + 9.*Fhydro.y*ux2*uy2 + 1/3.*Fhydro.y)/rho;\n' \
                          'F_in_cm[7] = m00*(9.*Fhydro.x*ux2*uy2 + 1/3.*Fhydro.x + 9.*Fhydro.y*u.x*uy3)/rho;\n' \
                          'F_in_cm[8] = -m00*(18.*Fhydro.x*ux3*uy2 + Fhydro.x*ux3 + 3.*Fhydro.x*u.x*uy2 + 18.*Fhydro.y*ux2*uy3 + 3.*Fhydro.y*ux2*u.y + Fhydro.y*uy3)/rho;\n'  # noqa

        assert 'F_in_cm[0] = 0;' in out
        assert 'F_in_cm[1] = Fhydro.x*m00/rho;' in out
        assert 'F_in_cm[2] = Fhydro.y*m00/rho;' in out
        assert 'F_in_cm[3] = -3.*m00*ux2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' in out
        assert 'F_in_cm[4] = -3.*m00*uy2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' in out
        assert 'F_in_cm[5] = -3.*m00*uxuy*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' in out
        assert 'F_in_cm[6] = m00*(9.*Fhydro.x*ux3*u.y + 9.*Fhydro.y*ux2*uy2 + 1/3.*Fhydro.y)/rho;\n' in out
        assert 'F_in_cm[7] = m00*(9.*Fhydro.x*ux2*uy2 + 1/3.*Fhydro.x + 9.*Fhydro.y*u.x*uy3)/rho;\n' in out
        assert 'F_in_cm[8] = -m00*(18.*Fhydro.x*ux3*uy2 + Fhydro.x*ux3 + 3.*Fhydro.x*u.x*uy2 + 18.*Fhydro.y*ux2*uy3 + 3.*Fhydro.y*ux2*u.y + Fhydro.y*uy3)/rho;\n' in out  # noqa

        assert expected_result == out

    def test_get_cm_eq_hydro_discrete(self):
        cm_eq = get_mom_vector_from_discrete_def(get_discrete_EDF_hydro, discrete_transform=get_discrete_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq', regex=True)
        out = f.getvalue()

        expected_result = 'cm_eq[0] = m00;\n' \
                          'cm_eq[1] = u.x*(-m00 + 1);\n' \
                          'cm_eq[2] = u.y*(-m00 + 1);\n' \
                          'cm_eq[3] = m00*ux2 + 1/3.*m00 - ux2;\n' \
                          'cm_eq[4] = m00*uy2 + 1/3.*m00 - uy2;\n' \
                          'cm_eq[5] = uxuy*(m00 - 1.);\n' \
                          'cm_eq[6] = u.y*(-m00*ux2 - 1/3.*m00 + 1/3.);\n' \
                          'cm_eq[7] = u.x*(-m00*uy2 - 1/3.*m00 + 1/3.);\n' \
                          'cm_eq[8] = m00*ux2*uy2 + 1/3.*m00*ux2 + 1/3.*m00*uy2 + 1/9.*m00 + 2.*ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n'  # noqa

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[1] = u.x*(-m00 + 1)' in out
        assert 'cm_eq[2] = u.y*(-m00 + 1);' in out
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

        cm_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * get_gamma(i), discrete_transform=get_discrete_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq', regex=True)
        out = f.getvalue()

        expected_result = 'cm_eq[0] = m00;\n' \
                          'cm_eq[1] = 0;\n' \
                          'cm_eq[2] = 0;\n' \
                          'cm_eq[3] = 1/3.*m00;\n' \
                          'cm_eq[4] = 1/3.*m00;\n' \
                          'cm_eq[5] = 0;\n' \
                          'cm_eq[6] = -m00*ux2*u.y;\n' \
                          'cm_eq[7] = -m00*u.x*uy2;\n' \
                          'cm_eq[8] = m00*(3.*ux2*uy2 + 1/9.);\n'

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
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSymbolicCalc)
    concurrent_suite = ConcurrentTestSuite(suite, fork_for_tests(cores))
    runner.run(concurrent_suite)
