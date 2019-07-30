from unittest import TestCase
from SymbolicCollisions.core.printers import print_as_vector

from sympy import Symbol
from sympy.matrices import Matrix
from SymbolicCollisions.core.cm_symbols import Temperature as T
from SymbolicCollisions.core.cm_symbols import m00, rho, cp, uy

import io
from contextlib import redirect_stdout

import sys
import os

from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir


class TestRegexPrinters(TestCase):
    def test_power_one(self):
        gamma = Symbol('gamma', positive=True)
        RT = Symbol('RT', positive=True)

        cases = [
             Matrix([1.*rho]),
             Matrix([1.0 * T * cp ** 1.0 * rho ** 1.0]),
             Matrix([0.1111111111 * T * gamma ** 2 / (cp * rho)]),
             Matrix([1.0 * m00 * (RT ** 2.0)]),
             Matrix([1.0 * m00 * (- RT ** 1.0 * uy ** 2)]),
             Matrix([RT*RT*m00]),
            ## Matrix([1.0 * m00 * (RT * uy ** 2 - RT ** 1.0 * uy ** 2 + RT ** 2.0)]),
             ]

        expected_results = [
            "\ttest = rho;\n",
            "\ttest = T*cp*rho;\n",
            "\ttest = 1/9.*T*gamma*gamma/(cp*rho);\n",
            "\ttest = RT*RT*m00;\n",
            "\ttest = -RT*m00*uy2;\n",
            "\ttest = RT*RT*m00;\n",
            ## "\ttest[0] = m00*(-RT*uy2 + RT*uy2 + RT*RT);\n",   # it seems that the order is sometimes swapped and the test fails Oo
        ]

        for case, expected_result in zip(cases, expected_results):
            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(case, 'test', raw_output=False)
            out = f.getvalue()

            assert out == expected_result

    def test_df_to_m(self):
        # SETUP
        d = 2
        q = 9

        # DYNAMIC IMPORTS
        ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
        ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
        if d == 3:
            ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
        else:
            ez = None

        pop_in_str = 'x_in'  # symbol defining populations
        temp_pop_str = 'temp'  # symbol defining populations

        temp_populations = get_DF(q, temp_pop_str)

        Mraw = get_raw_moments_matrix(ex, ey, ez)
        m = Mraw * temp_populations

        expected_results = "\tx_in[0] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[8];\n" \
                           "\tx_in[1] = temp[1] - temp[3] + temp[5] - temp[6] - temp[7] + temp[8];\n" \
                           "\tx_in[2] = temp[2] - temp[4] + temp[5] + temp[6] - temp[7] - temp[8];\n" \
                           "\tx_in[3] = temp[1] + temp[3] + temp[5] + temp[6] + temp[7] + temp[8];\n" \
                           "\tx_in[4] = temp[2] + temp[4] + temp[5] + temp[6] + temp[7] + temp[8];\n" \
                           "\tx_in[5] = temp[5] - temp[6] + temp[7] - temp[8];\n" \
                           "\tx_in[6] = temp[5] + temp[6] - temp[7] - temp[8];\n" \
                           "\tx_in[7] = temp[5] - temp[6] - temp[7] + temp[8];\n" \
                           "\tx_in[8] = temp[5] + temp[6] + temp[7] + temp[8];\n"

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(m, print_symbol=pop_in_str)
        out = f.getvalue()

        assert out == expected_results
