from unittest import TestCase
from SymbolicCollisions.core.printers import print_as_vector

from sympy import Symbol
from sympy.matrices import Matrix
from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import rho
from SymbolicCollisions.core.cm_symbols import Temperature as T
from SymbolicCollisions.core.cm_symbols import cp, ux, uy

import io
from contextlib import redirect_stdout

import sys
import os

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir


class TestRegexPrinters(TestCase):
    def test_power_one(self):
        gamma = Symbol('gamma', positive=True)
        RT = Symbol('RT', positive=True)

        cases = [Matrix([1.0 * T * cp ** 1.0 * rho ** 1.0]),
                 Matrix([0.1111111111 * T * gamma ** 2 / (cp * rho)]),
                 Matrix([1.0 * m00 * (RT * uy ** 2 - RT ** 1.0 * uy ** 2 + RT ** 2.0)])
                 ]

        expected_results = ["\ttest[0] = T*cp*rho;\n",
                            "\ttest[0] = 1/9.*T*gamma*gamma/(cp*rho);\n",
                            "\ttest[0] = m00*(RT*uy2 - RT*uy2 + RT*RT);\n"]

        for case, expected_result in zip(cases, expected_results):
            f = io.StringIO()
            with redirect_stdout(f):
                print_as_vector(case, 'test', raw_output=False)
            out = f.getvalue()

            assert out == expected_result
