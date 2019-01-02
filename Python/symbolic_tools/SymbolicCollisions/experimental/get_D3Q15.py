
from SymbolicCollisions.core.printers import print_as_vector
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix
from SymbolicCollisions.core.cm_symbols import ex_D2Q9, ey_D2Q9, Mraw_D2Q9

# from SymbolicCollisions.core.cm_symbols import \
#     ex_D3Q7 as ex, \
#     ey_D3Q7 as ey, \
#     ez_D3Q7 as ez

from SymbolicCollisions.core.cm_symbols import \
    ex_D3Q15 as ex, \
    ey_D3Q15 as ey, \
    ez_D3Q15 as ez

# from SymbolicCollisions.core.cm_symbols import \
#     ex_D3Q19 as ex, \
#     ey_D3Q19 as ey, \
#     ez_D3Q19 as ez

# from SymbolicCollisions.core.cm_symbols import \
#     ex_D2Q9 as ex, \
#     ey_D2Q9 as ey
#
# ez = None

from sympy import diff, ln, sin, pprint

import sys, os
sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir


# M = MatrixGenerator().get_raw_moments_matrix(ex_=ex, ey_=ey, ez_=ez)
M = get_raw_moments_matrix(ex_=ex, ey_=ey, ez_=ez)
# print_as_vector(Mraw, 's', regex=True)
pprint(M)

Smat = get_shift_matrix(M.inv(), ex, ey, ez_=ez)
pprint(Smat)
print_as_vector(Smat, 'N', regex=True)

