from SymbolicCollisions.core.printers import print_as_vector
from sympy import Symbol
from sympy.matrices import Matrix
from SymbolicCollisions.core.printers import round_and_simplify
from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import rho
from SymbolicCollisions.core.cm_symbols import Temperature as T
from SymbolicCollisions.core.cm_symbols import cp, ux, uy

"""
Here are samples to facilitate debugging regex-printers ;)
"""

gamma = Symbol('gamma', positive=True)
RT = Symbol('RT', positive=True)
# example= Matrix([1.0 * T * cp ** 1.0 * rho ** 1.0])
# example = Matrix([0.1111111111 * T * gamma ** 2 / (cp * rho)])
# example = Matrix([1.0 * m00 * (RT * uy ** 2 - RT ** 1.0 * uy ** 2 + RT ** 2.0)])

example = Matrix([1.0 * m00 * uy * (RT ** 1.0 - RT)])
print_as_vector(example, 'abc', raw_output=False)
