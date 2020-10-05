import numpy as np
from unittest import TestCase
from sympy.matrices import Matrix
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, get_reverse_indices


import sys
import os

sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir

class TestNotation(TestCase):

    def test_notation_i(self):
        clip_z_dimension = True

        m_seed = [0, 1, 2]
        rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)

        e_seed = [0, 1, -1]
        ex, ey, ez, e_new = get_e_as_in_r(e_seed, e_seed, e_seed)

        if clip_z_dimension:
            rmoments_order = rmoments_order[0:9]
            q, d = rmoments_order.shape
            d = 2
            ex = ex[0:9]
            ey = ey[0:9]
            ez = ez[0:9]
            e_D2Q9 = e_new[0:9, :]
        else:
            q, d = rmoments_order.shape

        expected_e = Matrix(
                    [[0,  0, 0],
                     [1,  0, 0],
                     [-1, 0, 0],
                     [0,  1, 0],
                     [1,  1, 0],
                     [-1, 1, 0],
                     [0, -1, 0],
                     [1, -1, 0],
                     [-1, -1, 0]])

        assert e_D2Q9 == expected_e

        rev_i = get_reverse_indices(e_D2Q9)
        expected_rev_i = np.array([0, 2, 1, 6, 8, 7, 3, 5, 4])

        np.testing.assert_array_equal(rev_i, expected_rev_i)
