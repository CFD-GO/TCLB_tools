from unittest import TestCase
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_indx_notation

from sympy import Symbol
from sympy.matrices import Matrix
from SymbolicCollisions.core.cm_symbols import Temperature as T
from SymbolicCollisions.core.cm_symbols import m00, rho, cp, uy
from SymbolicCollisions.core.cm_symbols import e_D3Q27, moments_dict


import io
from contextlib import redirect_stdout

import sys
import os

from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir


class TestRegexPrinters(TestCase):

    def test_notation_i(self):
        q = 27
        populations = get_print_symbols_in_indx_notation(q, print_symbol='b')
        expected_results = "\ta[0] = b[0];\n" \
                           "\ta[1] = b[1];\n" \
                           "\ta[2] = b[2];\n" \
                           "\ta[3] = b[3];\n" \
                           "\ta[4] = b[4];\n" \
                           "\ta[5] = b[5];\n" \
                           "\ta[6] = b[6];\n" \
                           "\ta[7] = b[7];\n" \
                           "\ta[8] = b[8];\n" \
                           "\ta[9] = b[9];\n" \
                           "\ta[10] = b[10];\n" \
                           "\ta[11] = b[11];\n" \
                           "\ta[12] = b[12];\n" \
                           "\ta[13] = b[13];\n" \
                           "\ta[14] = b[14];\n" \
                           "\ta[15] = b[15];\n" \
                           "\ta[16] = b[16];\n" \
                           "\ta[17] = b[17];\n" \
                           "\ta[18] = b[18];\n" \
                           "\ta[19] = b[19];\n" \
                           "\ta[20] = b[20];\n" \
                           "\ta[21] = b[21];\n" \
                           "\ta[22] = b[22];\n" \
                           "\ta[23] = b[23];\n" \
                           "\ta[24] = b[24];\n" \
                           "\ta[25] = b[25];\n" \
                           "\ta[26] = b[26];\n" \

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(populations, outprint_symbol='a')
        out = f.getvalue()

        assert out == expected_results

    def test_notation_ijk(self):
        q = 27
        populations = get_print_symbols_in_indx_notation(q, print_symbol='b')
        expected_results = "\ta000 = b[0];\n" \
                           "\ta100 = b[1];\n" \
                           "\ta200 = b[2];\n" \
                           "\ta010 = b[3];\n" \
                           "\ta020 = b[4];\n" \
                           "\ta001 = b[5];\n" \
                           "\ta002 = b[6];\n" \
                           "\ta111 = b[7];\n" \
                           "\ta211 = b[8];\n" \
                           "\ta121 = b[9];\n" \
                           "\ta221 = b[10];\n" \
                           "\ta112 = b[11];\n" \
                           "\ta212 = b[12];\n" \
                           "\ta122 = b[13];\n" \
                           "\ta222 = b[14];\n" \
                           "\ta110 = b[15];\n" \
                           "\ta210 = b[16];\n" \
                           "\ta120 = b[17];\n" \
                           "\ta220 = b[18];\n" \
                           "\ta101 = b[19];\n" \
                           "\ta201 = b[20];\n" \
                           "\ta102 = b[21];\n" \
                           "\ta202 = b[22];\n" \
                           "\ta011 = b[23];\n" \
                           "\ta021 = b[24];\n" \
                           "\ta012 = b[25];\n" \
                           "\ta022 = b[26];\n"

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(populations, outprint_symbol='a', output_order_of_moments=e_D3Q27)
        out = f.getvalue()

        assert out == expected_results

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

        temp_populations = get_print_symbols_in_indx_notation(q, temp_pop_str)

        matrixGenerator = MatrixGenerator(ex, ey, ez, moments_dict[f'D{d}Q{q}'])
        Mraw = matrixGenerator.get_raw_moments_matrix()

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
            print_as_vector(m, outprint_symbol=pop_in_str)
        out = f.getvalue()

        assert out == expected_results
