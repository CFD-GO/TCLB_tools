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

from SymbolicCollisions.core.printers import print_as_vector

from SymbolicCollisions.core.hardcoded_results import hardcoded_F_cm_hydro_compressible_D3Q19, \
    hardcoded_F_cm_Guo_hydro_incompressible_D2Q9, \
    hardcoded_F_cm_hydro_compressible_D2Q9, \
    hardcoded_cm_eq_compressible_D2Q9,    hardcoded_cm_eq_compressible_D3Q19, \
    hardcoded_cm_eq_incompressible_D2Q9, \
    hardcoded_cm_eq_compressible_D2Q9_thermal, \
    hardcoded_cm_eq_cht_D2Q9

from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, \
    rho, cs2_thermal, \
    m00, Force_str

from SymbolicCollisions.core.cm_symbols import moments_dict


class TestContinousCMTransforms(unittest.TestCase):
    def test_get_cm_eq_incompressible_continuous(self):
        # population_eq -> cm_eq - from continous definition: '
        # k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) '
        # where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n * (z-uz)^o ')

        cm_i = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
        cm_eq = get_mom_vector_from_continuous_def(cm_i.get_incompressible_DF,
                                                   continuous_transformation=cm_i.get_cm,
                                                   moments_order=moments_dict['D2Q9'],
                                                   serial_run=True)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq')
        out = f.getvalue()

        # TODO: can't use hardcoded_cm_eq_incompressible_D2Q9,
        #  because sympy switches hardcoded 'u.x*(-m00 + 1)' to '-u.x*(m00 - 1') and test fails.
        #  thank you sympy...

        assert f'cm_eq[0] = {m00};' in out
        assert f'cm_eq[1] = u.x*(1 - {m00});' in out or f'cm_eq[1] = u.x*(-{m00} + 1);' in out
        assert f'cm_eq[2] = u.y*(1 - {m00});' in out or f'cm_eq[2] = u.y*(-{m00} + 1);' in out
        assert f'cm_eq[3] = {m00}*ux2 + 1/3.*{m00} - ux2;\n' in out
        assert f'cm_eq[4] = {m00}*uy2 + 1/3.*{m00} - uy2;\n' in out
        assert f'cm_eq[5] = uxuy*({m00} - 1.);\n' in out
        assert f'cm_eq[6] = u.y*(-{m00}*ux2 - 1/3.*{m00} + ux2 + 1/3.);\n' in out
        assert f'cm_eq[7] = u.x*(-{m00}*uy2 - 1/3.*{m00} + uy2 + 1/3.);\n' in out
        assert f'cm_eq[8] = {m00}*ux2*uy2 + 1/3.*{m00}*ux2 + 1/3.*{m00}*uy2 + 1/9.*{m00} - ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n' in out  # noqa


        # expected_result = f'\tcm_eq[0] = {m00};\n' \
        #                   f'\tcm_eq[1] = u.x*(1 - {m00});\n' \
        #                   f'\tcm_eq[2] = u.y*(1 - {m00});\n' \
        #                   f'\tcm_eq[3] = {m00}*ux2 + 1/3.*{m00} - ux2;\n' \
        #                   f'\tcm_eq[4] = {m00}*uy2 + 1/3.*{m00} - uy2;\n' \
        #                   f'\tcm_eq[5] = uxuy*({m00} - 1.);\n' \
        #                   f'\tcm_eq[6] = u.y*(-{m00}*ux2 - 1/3.*{m00} + ux2 + 1/3.);\n' \
        #                   f'\tcm_eq[7] = u.x*(-{m00}*uy2 - 1/3.*{m00} + uy2 + 1/3.);\n' \
        #                   f'\tcm_eq[8] = {m00}*ux2*uy2 + 1/3.*{m00}*ux2 + 1/3.*{m00}*uy2 + 1/9.*{m00} - ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n'  # noqa

        # assert expected_result == out

    def test_get_F_cm_using_He_scheme_and_continuous_Maxwellian_DF(self):
        cm_i = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
        F_cm = get_mom_vector_from_continuous_def(cm_i.get_force_He_hydro_DF,
                                                  continuous_transformation=cm_i.get_cm,
                                                  moments_order=moments_dict['D2Q9'],
                                                  serial_run=True)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_cm, 'F_cm')
        out = f.getvalue()

        # TODO: can't use hardcoded_cm_eq_incompressible_D2Q9,
        #  because sympy switches hardcoded terms like 'u.x*(-m00 + 1)' to '-u.x*(m00 - 1') and test fails.
        #  thank you sympy...


        assert f'\tF_cm[0] = 0;\n' in out
        assert f'\tF_cm[1] = {Force_str}.x*{m00}/rho;\n' in out
        assert f'\tF_cm[2] = {Force_str}.y*{m00}/rho;\n' in out
        assert f'\tF_cm[3] = -2.*{Force_str}.x*u.x*({m00} - 1.)/rho;\n' in out
        assert f'\tF_cm[4] = -2.*{Force_str}.y*u.y*({m00} - 1.)/rho;\n' in out
        assert f'\tF_cm[5] = (-{Force_str}.x*{m00}*u.y + {Force_str}.x*u.y - {Force_str}.y*{m00}*u.x + {Force_str}.y*u.x)/rho;\n' in out
        assert f'\tF_cm[6] = (2.*{Force_str}.x*{m00}*uxuy - 2.*{Force_str}.x*uxuy + {Force_str}.y*{m00}*ux2 + 1/3.*{Force_str}.y*{m00} - {Force_str}.y*ux2)/rho;\n' in out
        assert f'\tF_cm[7] = ({Force_str}.x*{m00}*uy2 + 1/3.*{Force_str}.x*{m00} - {Force_str}.x*uy2 + 2.*{Force_str}.y*{m00}*uxuy - 2.*{Force_str}.y*uxuy)/rho;\n' in out
        assert f'\tF_cm[8] = (-2.*{Force_str}.x*{m00}*u.x*uy2 - 2/3.*{Force_str}.x*{m00}*u.x + 2.*{Force_str}.x*u.x*uy2 + 2/3.*{Force_str}.x*u.x - 2.*{Force_str}.y*{m00}*ux2*u.y - 2/3.*{Force_str}.y*{m00}*u.y + 2.*{Force_str}.y*ux2*u.y + 2/3.*{Force_str}.y*u.y)/rho;\n' in out


        # expected_result = \
        #     f'\tF_cm[0] = 0;\n' \
        #     f'\tF_cm[1] = {Force_str}.x*{m00}/rho;\n' \
        #     f'\tF_cm[2] = {Force_str}.y*{m00}/rho;\n' \
        #     f'\tF_cm[3] = -2.*{Force_str}.x*u.x*({m00} - 1.)/rho;\n' \
        #     f'\tF_cm[4] = -2.*{Force_str}.y*u.y*({m00} - 1.)/rho;\n' \
        #     f'\tF_cm[5] = (-{Force_str}.x*{m00}*u.y + {Force_str}.x*u.y - {Force_str}.y*{m00}*u.x + {Force_str}.y*u.x)/rho;\n' \
        #     f'\tF_cm[6] = (2.*{Force_str}.x*{m00}*uxuy - 2.*{Force_str}.x*uxuy + {Force_str}.y*{m00}*ux2 + 1/3.*{Force_str}.y*{m00} - {Force_str}.y*ux2)/rho;\n' \
        #     f'\tF_cm[7] = ({Force_str}.x*{m00}*uy2 + 1/3.*{Force_str}.x*{m00} - {Force_str}.x*uy2 + 2.*{Force_str}.y*{m00}*uxuy - 2.*{Force_str}.y*uxuy)/rho;\n' \
        #     f'\tF_cm[8] = (-2.*{Force_str}.x*{m00}*u.x*uy2 - 2/3.*{Force_str}.x*m00*u.x + 2.*{Force_str}.x*u.x*uy2 + 2/3.*{Force_str}.x*u.x - 2.*{Force_str}.y*{m00}*ux2*u.y - 2/3.*{Force_str}.y*{m00}*u.y + 2.*{Force_str}.y*ux2*u.y + 2/3.*{Force_str}.y*u.y)/rho;\n'
        #
        # assert expected_result == out

    def test_thermal_cm_eq_vector_from_continuous_def(self):
        ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho, cs2=cs2_thermal)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                                      continuous_transformation=ccmt.get_cm,
                                                      moments_order=moments_dict['D2Q9'],
                                                      serial_run=True)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq')
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            expected_result = hardcoded_cm_eq_compressible_D2Q9_thermal
            print_as_vector(expected_result, 'cm_eq')

        ccode_expected_result = f2.getvalue()

        assert ccode_expected_result == out

    def test_cm_vector_from_continuous_def(self):
        # this test runs long without output and CI may consider it as a timeout :/
        ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)

        from SymbolicCollisions.core.cm_symbols import Sigma2asSymbol
        ccmt_cht = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho, cs2=Sigma2asSymbol)

        ccmts = [
            ccmt,
            ccmt,
            ccmt,
            ccmt,
            ccmt,
            ccmt_cht,
        ]

        lattices = [
            'D2Q9',
            'D3Q19',
            'D2Q9',
            'D2Q9',
            'D3Q19',
            'D2Q9'
            ]

        functions = [
            ccmt.get_Maxwellian_DF,
            ccmt.get_Maxwellian_DF,
            ccmt.get_force_Guo,
            ccmt.get_force_He_MB,
            ccmt.get_force_He_MB,
            ccmt_cht.get_cht_DF,
        ]

        expected_results = [
            hardcoded_cm_eq_compressible_D2Q9,
            hardcoded_cm_eq_compressible_D3Q19,
            hardcoded_F_cm_Guo_hydro_incompressible_D2Q9,
            hardcoded_F_cm_hydro_compressible_D2Q9,
            hardcoded_F_cm_hydro_compressible_D3Q19,
            hardcoded_cm_eq_cht_D2Q9,
        ]

        for fun, lattice, _ccmt, expected_result in zip(functions, lattices, ccmts, expected_results):
            cm_eq = get_mom_vector_from_continuous_def(fun,
                                                       continuous_transformation=_ccmt.get_cm,
                                                       moments_order=moments_dict[lattice],
                                                       serial_run=True
                                                       )
            # print("------------\n\n")
            # print_as_vector(cm_eq, 'CM')
            # print_as_vector(expected_result, 'CM_expected')
            # print("------------\n\n")

            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(cm_eq, 'cm_eq')
            out = f.getvalue()

            f2 = io.StringIO()
            with redirect_stdout(f2):
                print_as_vector(expected_result, 'cm_eq')
            ccode_expected_result = f2.getvalue()

            assert ccode_expected_result == out


# Pycharm runs them sequentially
# python -m unittest tests/test_example_unit_tests_parallel_run.py # sequential as well
# python tests/test_example_unit_tests_parallel_run.py # concurrent run :)

if __name__ == '__main__':
    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner()
    # Run same tests across 4 processes
    cores = multiprocessing.cpu_count()
    print(f'\nRunning tests on {cores} cores:')
    suite = unittest.TestLoader().loadTestsFromTestCase(TestContinousCMTransforms)
    concurrent_suite = ConcurrentTestSuite(suite, fork_for_tests(cores))
    runner.run(concurrent_suite)
