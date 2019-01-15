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

from SymbolicCollisions.core.hardcoded_results import hardcoded_F_cm_hydro_density_based_D3Q19, \
    hardcoded_F_cm_Guo_hydro_LB_velocity_based_D2Q9, \
    hardcoded_F_cm_hydro_density_based_D2Q9, \
    hardcoded_cm_eq_compressible_D2Q9,    hardcoded_cm_eq_compressible_D3Q19, \
    hardcoded_cm_eq_incompressible_D2Q9

from SymbolicCollisions.core.ContinousCMTransforms import ContinousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, \
    F2D, dzeta2D, u2D, \
    rho, w_D2Q9, m00

from SymbolicCollisions.core.cm_symbols import moments_dict


class TestSymbolicCalc(unittest.TestCase):
    def test_cm_vector_from_continuous_def(self):
        ccmt = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)

        lattices = [
            'D2Q9',
            'D3Q19',
            'D2Q9',
            'D2Q9',
            'D3Q19',
            ]

        functions = [
            ccmt.get_Maxwellian_DF,
            ccmt.get_Maxwellian_DF,
            ccmt.get_force_Guo,
            ccmt.get_force_He_MB,
            ccmt.get_force_He_MB,
        ]

        expected_results = [
            hardcoded_cm_eq_compressible_D2Q9,
            hardcoded_cm_eq_compressible_D3Q19,
            hardcoded_F_cm_Guo_hydro_LB_velocity_based_D2Q9,
            hardcoded_F_cm_hydro_density_based_D2Q9,
            hardcoded_F_cm_hydro_density_based_D3Q19,
        ]

        for fun, lattice, expected_result in zip(functions, lattices, expected_results):
            cm_eq = get_mom_vector_from_continuous_def(fun,
                                                       continuous_transformation=ccmt.get_cm,
                                                       moments_order=moments_dict[lattice])
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

    def test_get_cm_eq_incompressible_continuous(self):
        # population_eq -> cm_eq - from continous definition: '
        # k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) '
        # where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n * (z-uz)^o ')

        cm_i = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
        cm_eq = get_mom_vector_from_continuous_def(cm_i.get_hydro_DF,
                                                   continuous_transformation=cm_i.get_cm,
                                                   moments_order=moments_dict['D2Q9'])

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(cm_eq, 'cm_eq')
        out = f.getvalue()

        # TODO: can't use hardcoded_cm_eq_incompressible_D2Q9,
        #  because sympy switches hardcoded 'u.x*(-m00 + 1)' to '-u.x*(m00 - 1') and test fails.
        #  thank you sympy...

        expected_result = 'cm_eq[0] = m00;\n' \
                          'cm_eq[1] = u.x*(-m00 + 1);\n' \
                          'cm_eq[2] = u.y*(-m00 + 1);\n' \
                          'cm_eq[3] = m00*ux2 + 1/3.*m00 - ux2;\n' \
                          'cm_eq[4] = m00*uy2 + 1/3.*m00 - uy2;\n' \
                          'cm_eq[5] = uxuy*(m00 - 1.);\n' \
                          'cm_eq[6] = u.y*(-m00*ux2 - 1/3.*m00 + ux2 + 1/3.);\n' \
                          'cm_eq[7] = u.x*(-m00*uy2 - 1/3.*m00 + uy2 + 1/3.);\n' \
                          'cm_eq[8] = m00*ux2*uy2 + 1/3.*m00*ux2 + 1/3.*m00*uy2 + 1/9.*m00 - ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n'  # noqa

        assert 'cm_eq[0] = m00;' in out
        assert 'cm_eq[1] = u.x*(-m00 + 1)' in out
        assert 'cm_eq[2] = u.y*(-m00 + 1);' in out
        assert 'cm_eq[3] = m00*ux2 + 1/3.*m00 - ux2;\n' in out
        assert 'cm_eq[4] = m00*uy2 + 1/3.*m00 - uy2;\n' in out
        assert 'cm_eq[5] = uxuy*(m00 - 1.);\n' in out
        assert 'cm_eq[6] = u.y*(-m00*ux2 - 1/3.*m00 + ux2 + 1/3.);\n' in out
        assert 'cm_eq[7] = u.x*(-m00*uy2 - 1/3.*m00 + uy2 + 1/3.);\n' in out
        assert 'cm_eq[8] = m00*ux2*uy2 + 1/3.*m00*ux2 + 1/3.*m00*uy2 + 1/9.*m00 - ux2*uy2 - 1/3.*ux2 - 1/3.*uy2;\n' in out  # noqa

        assert expected_result == out

    def test_get_F_cm_using_He_scheme_and_continuous_Maxwellian_DF(self):
        cm_i = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
        F_cm = get_mom_vector_from_continuous_def(cm_i.get_force_He_hydro_DF,
                                                  continuous_transformation=cm_i.get_cm,
                                                  moments_order=moments_dict['D2Q9'])

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(F_cm, 'F_cm')
        out = f.getvalue()

        # TODO: can't use hardcoded_cm_eq_incompressible_D2Q9,
        #  because sympy switches hardcoded terms like 'u.x*(-m00 + 1)' to '-u.x*(m00 - 1') and test fails.
        #  thank you sympy...

        expected_result = \
            'F_cm[0] = 0;\n' \
            'F_cm[1] = Fhydro.x*m00/rho;\n' \
            'F_cm[2] = Fhydro.y*m00/rho;\n' \
            'F_cm[3] = -2.*Fhydro.x*u.x*(m00 - 1.)/rho;\n' \
            'F_cm[4] = -2.*Fhydro.y*u.y*(m00 - 1.)/rho;\n' \
            'F_cm[5] = (-Fhydro.x*m00*u.y + Fhydro.x*u.y - Fhydro.y*m00*u.x + Fhydro.y*u.x)/rho;\n' \
            'F_cm[6] = (2.*Fhydro.x*m00*uxuy - 2.*Fhydro.x*uxuy + Fhydro.y*m00*ux2 + 1/3.*Fhydro.y*m00 - Fhydro.y*ux2)/rho;\n' \
            'F_cm[7] = (Fhydro.x*m00*uy2 + 1/3.*Fhydro.x*m00 - Fhydro.x*uy2 + 2.*Fhydro.y*m00*uxuy - 2.*Fhydro.y*uxuy)/rho;\n' \
            'F_cm[8] = (-2.*Fhydro.x*m00*u.x*uy2 - 2/3.*Fhydro.x*m00*u.x + 2.*Fhydro.x*u.x*uy2 + 2/3.*Fhydro.x*u.x - 2.*Fhydro.y*m00*ux2*u.y - 2/3.*Fhydro.y*m00*u.y + 2.*Fhydro.y*ux2*u.y + 2/3.*Fhydro.y*u.y)/rho;\n'

        assert expected_result == out
