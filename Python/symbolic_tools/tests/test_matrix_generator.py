from unittest import TestCase

import io
from contextlib import redirect_stdout

import sys
import os
sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir

import unittest
from SymbolicCollisions.core.printers import print_as_vector

import multiprocessing
from concurrencytest import ConcurrentTestSuite, fork_for_tests
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator
from SymbolicCollisions.core.cm_symbols import moments_dict

class TestMatrixGenerator(TestCase):

    def test_get_raw_matrix_d2q9(self):

        from SymbolicCollisions.core.cm_symbols import ex_D2Q9, ey_D2Q9, Mraw_D2Q9

        matrix_generator = MatrixGenerator(ex=ex_D2Q9, ey=ey_D2Q9, ez=None, order_of_moments=moments_dict['D2Q9'])
        M = matrix_generator.get_raw_moments_matrix()

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(Mraw_D2Q9, 's')
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(M, 's')
        out2 = f2.getvalue()

        assert out == out2

    def test_Shift_ortho_Straka_d2q5(self):
        from SymbolicCollisions.core.cm_symbols import Shift_ortho_Straka_d2q5, K_ortho_Straka_d2q5, \
            ex_Straka_d2_q5, ey_Straka_d2_q5

        matrix_generator = MatrixGenerator(ex=ex_Straka_d2_q5, ey=ey_Straka_d2_q5, ez=None, order_of_moments=moments_dict['D2Q5'])
        Smat = matrix_generator.get_shift_matrix(K_ortho_Straka_d2q5)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(Shift_ortho_Straka_d2q5, 's')
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(Smat[1:, 1:], 's')
        out2 = f2.getvalue()

        assert out == out2

    def test_Shift_ortho_Geier_d2q9(self):

        from SymbolicCollisions.core.cm_symbols import Shift_ortho_Geier, K_ortho_Geier, ex_Geier, ey_Geier
        matrix_generator = MatrixGenerator(ex=ex_Geier, ey=ey_Geier, ez=None, order_of_moments=moments_dict['D2Q9'])
        Smat = matrix_generator.get_shift_matrix(K_ortho_Geier)

        f = io.StringIO()
        with redirect_stdout(f):
            print_as_vector(Shift_ortho_Geier, 's')
        out = f.getvalue()

        f2 = io.StringIO()
        with redirect_stdout(f2):
            print_as_vector(Smat[3:, 3:], 's')
        out2 = f2.getvalue()

        assert out == out2

# Pycharm runs them sequentially
# python -m unittest tests/test_example_unit_tests_parallel_run.py # sequential as well
# python tests/test_example_unit_tests_parallel_run.py # concurrent run :)

if __name__ == '__main__':
    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner()
    # Run same tests across 4 processes
    cores = multiprocessing.cpu_count()
    print(f'\nRunning tests on {cores} cores:')
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMatrixGenerator)
    concurrent_suite = ConcurrentTestSuite(suite, fork_for_tests(cores))
    # concurrent_suite.addTests(loader.loadTestsFromTestCase(TestMatrixGenerator))  # add more :>
    runner.run(concurrent_suite)
