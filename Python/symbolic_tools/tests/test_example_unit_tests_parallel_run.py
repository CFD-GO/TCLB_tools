import time
import unittest

import multiprocessing
from concurrencytest import ConcurrentTestSuite, fork_for_tests


class SampleTestCase(unittest.TestCase):
    """Dummy tests that sleep for demo."""

    def test_me_1(self):
        time.sleep(0.25)

    def test_me_2(self):
        time.sleep(0.25)

    def test_me_3(self):
        time.sleep(0.25)
        # assert False

    def test_me_4(self):
        time.sleep(0.25)


# Pycharm runs them sequentially
# python -m unittest tests/test_example_unit_tests_parallel_run.py # sequential as well
# python tests/test_example_unit_tests_parallel_run.py # concurrent run :)

if __name__ == '__main__':
    # Load tests from SampleTestCase defined above
    loader = unittest.TestLoader()

    suite = loader.loadTestsFromTestCase(SampleTestCase)
    suite.addTests(loader.loadTestsFromTestCase(SampleTestCase))  # add more :>
    runner = unittest.TextTestRunner()

    # Run tests sequentially
    print('\nRun tests sequentially:')
    runner.run(suite)

    # Run same tests across 4 processes
    cores = multiprocessing.cpu_count()
    print(f'\nRun tests on {cores} cores:')
    suite = unittest.TestLoader().loadTestsFromTestCase(SampleTestCase)
    concurrent_suite = ConcurrentTestSuite(suite, fork_for_tests(cores))
    # concurrent_suite.addTests(loader.loadTestsFromTestCase(SampleTestCase))  # add more :>
    runner.run(concurrent_suite)

