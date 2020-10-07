from unittest import TestCase

import io
from contextlib import redirect_stdout

import sys
import os
sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir


class TestCumulants(TestCase):

    def test_get_cumulants(self):
        from SymbolicCollisions.core.cumulants import get_cumulant
        from SymbolicCollisions.core.printers import print_as_vector
        from sympy.matrices import Matrix

        order = [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (2, 0, 0),
            (0, 2, 0),
            (0, 0, 2),
            (1, 1, 0)]

        result = [
            '= m100;',
            '= m010;',
            '= m001;',
            '= m200 - m100*m100/m000;',
            '= m020 - m010*m010/m000',
            '= m002 - m001*m001/m000;',
            '= m110 - m010*m100/m000;']

        for o, r in zip(order, result):
            cumulant, _ = get_cumulant(*o)

            f = io.StringIO()
            with redirect_stdout(f):
                mc = Matrix([cumulant])  # printer works using Matrix format
                print_as_vector(mc, 'c')
            out = f.getvalue()

            assert r in out
