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

from SymbolicCollisions.core.sym_col_fun import \
    get_mom_vector_from_discrete_def, get_mom_vector_from_shift_Mat, \
    get_mom_vector_from_continuous_def, get_continuous_Maxwellian_DF, \
    get_continuous_force_He_MB, \
    get_discrete_EDF_hydro, \
    get_discrete_force_He, \
    get_discrete_force_Guo, get_continuous_force_Guo, \
    get_gamma, get_continuous_hydro_DF, get_continuous_force_He_hydro_DF, \
    get_continuous_cm, get_discrete_cm

from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, NrawD2Q9

from SymbolicCollisions.core.printers import print_as_vector

from SymbolicCollisions.core.hardcoded_results import \
    hardcoded_F_cm_Guo_hydro_LB_velocity_based, hardcoded_cm_pf_eq, hardcoded_cm_hydro_eq


class TestSymbolicCalc(unittest.TestCase):
    def test_get_raw_matrix_d2q9(self):
        from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix
        from SymbolicCollisions.core.cm_symbols import ex_D2Q9, ey_D2Q9, Mraw_D2Q9

        M = get_raw_moments_matrix(ex_=ex_D2Q9, ey_=ey_D2Q9)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(Mraw_D2Q9, 's', regex=True)
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(M, 's', regex=True)
        out2 = f2.getvalue()

        assert out == out2


    def test_Shift_ortho_Straka_d2q5(self):
        from SymbolicCollisions.core.MatrixGenerator import get_shift_matrix
        from SymbolicCollisions.core.cm_symbols import Shift_ortho_Straka_d2q5, K_ortho_Straka_d2q5, ex_Straka_d2_q5, \
            ey_Straka_d2_q5

        Smat = get_shift_matrix(K_ortho_Straka_d2q5, ex_Straka_d2_q5, ey_Straka_d2_q5)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(Shift_ortho_Straka_d2q5, 's', regex=True)
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(Smat[1:, 1:], 's', regex=True)
        out2 = f2.getvalue()

        assert out == out2

    def test_Shift_ortho_Geier_d2q9(self):
        from SymbolicCollisions.core.MatrixGenerator import get_shift_matrix
        from SymbolicCollisions.core.cm_symbols import Shift_ortho_Geier, K_ortho_Geier, ex_Geier, ey_Geier

        Smat = get_shift_matrix(K_ortho_Geier, ex_Geier, ey_Geier)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(Shift_ortho_Geier, 's', regex=True)
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(Smat[3:, 3:], 's', regex=True)
        out2 = f2.getvalue()

        assert out == out2

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

    def test_get_cm_eq_from_continuous_Maxwellian_DF(self):
        cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq', regex=True)
        out = f.getvalue()

        f2= io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(hardcoded_cm_pf_eq, 'cm_eq', regex=True)
        expected_result = f2.getvalue()

        assert expected_result == out

    def test_get_F_cm_Guo_continuous_and_discrete(self):
        F_cm_Guo_disc = get_mom_vector_from_discrete_def(get_discrete_force_Guo, discrete_transform=get_discrete_cm)
        F_cm_Guo_cont = get_mom_vector_from_continuous_def(get_continuous_force_Guo, continuous_transformation=get_continuous_cm)

        results = [F_cm_Guo_disc, F_cm_Guo_cont]

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(hardcoded_F_cm_Guo_hydro_LB_velocity_based, 'F_cm', regex=True)
        expected_result = f.getvalue()

        for result in results:
            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(result, 'F_cm', regex=True)
            out = f.getvalue()

        assert out == expected_result

    def test_get_F_cm_using_He_scheme_and_continuous_rho_Maxwellian_DF(self):
        from SymbolicCollisions.core.hardcoded_results import hardcoded_F_cm_hydro_LB_density_based

        F_cm = get_mom_vector_from_continuous_def(get_continuous_force_He_MB, continuous_transformation=get_continuous_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_cm, 'F_cm', regex=True)
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(hardcoded_F_cm_hydro_LB_density_based, 'F_cm', regex=True)
            expected_result = f2.getvalue()

        assert expected_result == out

    def test_get_F_cm_using_He_scheme_and_continuous_p_Maxwellian_DF(self):
        F_cm = get_mom_vector_from_continuous_def(get_continuous_force_He_hydro_DF, continuous_transformation=get_continuous_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_cm, 'F_cm', regex=True)
        out = f.getvalue()

        expected_result = 'F_cm[0] = 0;\n' \
                          'F_cm[1] = Fhydro.x*m00/rho;\n' \
                          'F_cm[2] = Fhydro.y*m00/rho;\n' \
                          'F_cm[3] = -2.0*Fhydro.x*u.x*(m00 - 1.0)/rho;\n' \
                          'F_cm[4] = -2.0*Fhydro.y*u.y*(m00 - 1.0)/rho;\n' \
                          'F_cm[5] = (-Fhydro.x*m00*u.y + Fhydro.x*u.y - Fhydro.y*m00*u.x + Fhydro.y*u.x)/rho;\n' \
                          'F_cm[6] = (2.0*Fhydro.x*m00*uxuy - 2.0*Fhydro.x*uxuy + Fhydro.y*m00*ux2 + 1./3.*Fhydro.y*m00 - Fhydro.y*ux2)/rho;\n' \
                          'F_cm[7] = (Fhydro.x*m00*uy2 + 1./3.*Fhydro.x*m00 - Fhydro.x*uy2 + 2.0*Fhydro.y*m00*uxuy - 2.0*Fhydro.y*uxuy)/rho;\n' \
                          'F_cm[8] = (-2.0*Fhydro.x*m00*u.x*uy2 - 2./3.*Fhydro.x*m00*u.x + 2.0*Fhydro.x*u.x*uy2 + 2./3.*Fhydro.x*u.x - 2.0*Fhydro.y*m00*ux2*u.y - 2./3.*Fhydro.y*m00*u.y + 2.0*Fhydro.y*ux2*u.y + 2./3.*Fhydro.y*u.y)/rho;\n'

        assert expected_result == out

    def test_get_force_He_original(self):
        F_in_cm = get_mom_vector_from_discrete_def(get_discrete_force_He, discrete_transform=get_discrete_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_in_cm, 'F_in_cm', regex=True)
        out = f.getvalue()

        expected_result = 'F_in_cm[0] = 0;\n' \
                          'F_in_cm[1] = Fhydro.x*m00/rho;\n' \
                          'F_in_cm[2] = Fhydro.y*m00/rho;\n' \
                          'F_in_cm[3] = -3.0*m00*ux2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' \
                          'F_in_cm[4] = -3.0*m00*uy2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' \
                          'F_in_cm[5] = -3.0*m00*uxuy*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' \
                          'F_in_cm[6] = m00*(9.0*Fhydro.x*ux3*u.y + 9.0*Fhydro.y*ux2*uy2 + 1./3.*Fhydro.y)/rho;\n' \
                          'F_in_cm[7] = m00*(9.0*Fhydro.x*ux2*uy2 + 1./3.*Fhydro.x + 9.0*Fhydro.y*uxuy3)/rho;\n' \
                          'F_in_cm[8] = -m00*(18.0*Fhydro.x*ux3*uy2 + Fhydro.x*ux3 + 3.0*Fhydro.x*u.x*uy2 + 18.0*Fhydro.y*ux2*uy3 + 3.0*Fhydro.y*ux2*u.y + Fhydro.y*uy3)/rho;\n'  # noqa

        assert 'F_in_cm[0] = 0;' in out
        assert 'F_in_cm[1] = Fhydro.x*m00/rho;' in out
        assert 'F_in_cm[2] = Fhydro.y*m00/rho;' in out
        assert 'F_in_cm[3] = -3.0*m00*ux2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' in out
        assert 'F_in_cm[4] = -3.0*m00*uy2*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' in out
        assert 'F_in_cm[5] = -3.0*m00*uxuy*(Fhydro.x*u.x + Fhydro.y*u.y)/rho;\n' in out
        assert 'F_in_cm[6] = m00*(9.0*Fhydro.x*ux3*u.y + 9.0*Fhydro.y*ux2*uy2 + 1./3.*Fhydro.y)/rho;\n' in out
        assert 'F_in_cm[7] = m00*(9.0*Fhydro.x*ux2*uy2 + 1./3.*Fhydro.x + 9.0*Fhydro.y*uxuy3)/rho;\n' in out
        assert 'F_in_cm[8] = -m00*(18.0*Fhydro.x*ux3*uy2 + Fhydro.x*ux3 + 3.0*Fhydro.x*u.x*uy2 + 18.0*Fhydro.y*ux2*uy3 + 3.0*Fhydro.y*ux2*u.y + Fhydro.y*uy3)/rho;\n' in out  # noqa

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
                          'cm_eq[3] = m00*ux2 + 1./3.*m00 - ux2;\n' \
                          'cm_eq[4] = m00*uy2 + 1./3.*m00 - uy2;\n' \
                          'cm_eq[5] = uxuy*(m00 - 1.0);\n' \
                          'cm_eq[6] = u.y*(-m00*ux2 - 1./3.*m00 + 1./3.);\n' \
                          'cm_eq[7] = u.x*(-m00*uy2 - 1./3.*m00 + 1./3.);\n' \
                          'cm_eq[8] = m00*ux2*uy2 + 1./3.*m00*ux2 + 1./3.*m00*uy2 + 1./9.*m00 + 2.0*ux2*uy2 - 1./3.*ux2 - 1./3.*uy2;\n'  # noqa

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[1] = u.x*(-m00 + 1)' in out
        assert 'cm_eq[2] = u.y*(-m00 + 1);' in out
        assert 'cm_eq[3] = m00*ux2 + 1./3.*m00 - ux2;\n' in out
        assert 'cm_eq[4] = m00*uy2 + 1./3.*m00 - uy2;\n' in out
        assert 'cm_eq[5] = uxuy*(m00 - 1.0);\n' in out
        assert 'cm_eq[6] = u.y*(-m00*ux2 - 1./3.*m00 + 1./3.);\n' in out
        assert 'cm_eq[7] = u.x*(-m00*uy2 - 1./3.*m00 + 1./3.);\n' in out
        assert 'cm_eq[8] = m00*ux2*uy2 + 1./3.*m00*ux2 + 1./3.*m00*uy2 + 1./9.*m00 + 2.0*ux2*uy2 - 1./3.*ux2 - 1./3.*uy2;\n' in out  # noqa

        assert expected_result == out

    def test_get_cm_eq_hydro_cont(self):
        # population_eq -> cm_eq - from continous definition: '
        # k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) '
        # where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')

        cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq', regex=True)
        out = f.getvalue()

        expected_result = 'cm_eq[0] = m00;\n' \
                          'cm_eq[1] = u.x*(-m00 + 1);\n' \
                          'cm_eq[2] = u.y*(-m00 + 1);\n' \
                          'cm_eq[3] = m00*ux2 + 1./3.*m00 - ux2;\n' \
                          'cm_eq[4] = m00*uy2 + 1./3.*m00 - uy2;\n' \
                          'cm_eq[5] = uxuy*(m00 - 1.0);\n' \
                          'cm_eq[6] = u.y*(-m00*ux2 - 1./3.*m00 + ux2 + 1./3.);\n' \
                          'cm_eq[7] = u.x*(-m00*uy2 - 1./3.*m00 + uy2 + 1./3.);\n' \
                          'cm_eq[8] = m00*ux2*uy2 + 1./3.*m00*ux2 + 1./3.*m00*uy2 + 1./9.*m00 - ux2*uy2 - 1./3.*ux2 - 1./3.*uy2;\n'  # noqa

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[1] = u.x*(-m00 + 1)' in out
        assert 'cm_eq[2] = u.y*(-m00 + 1);' in out
        assert 'cm_eq[3] = m00*ux2 + 1./3.*m00 - ux2;\n' in out
        assert 'cm_eq[4] = m00*uy2 + 1./3.*m00 - uy2;\n' in out
        assert 'cm_eq[5] = uxuy*(m00 - 1.0);\n' in out
        assert 'cm_eq[6] = u.y*(-m00*ux2 - 1./3.*m00 + ux2 + 1./3.);\n' in out
        assert 'cm_eq[7] = u.x*(-m00*uy2 - 1./3.*m00 + uy2 + 1./3.);\n' in out
        assert 'cm_eq[8] = m00*ux2*uy2 + 1./3.*m00*ux2 + 1./3.*m00*uy2 + 1./9.*m00 - ux2*uy2 - 1./3.*ux2 - 1./3.*uy2;\n' in out  # noqa

        assert expected_result == out

    def test_cm_eq(self):
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
                          'cm_eq[3] = 1./3.*m00;\n' \
                          'cm_eq[4] = 1./3.*m00;\n' \
                          'cm_eq[5] = 0;\n' \
                          'cm_eq[6] = -m00*ux2*u.y;\n' \
                          'cm_eq[7] = -m00*u.x*uy2;\n' \
                          'cm_eq[8] = m00*(3.0*ux2*uy2 + 1./9.);\n'

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[2] = 0;' in out
        assert 'cm_eq[2] = 0;' in out
        assert 'cm_eq[3] = 1./3.*m00;\n' in out
        assert 'cm_eq[4] = 1./3.*m00;\n' in out
        assert 'cm_eq[5] = 0;\n' in out
        assert 'cm_eq[6] = -m00*ux2*u.y;\n' in out
        assert 'cm_eq[7] = -m00*u.x*uy2;\n' in out
        assert 'cm_eq[8] = m00*(3.0*ux2*uy2 + 1./9.);\n' in out

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
