from sympy import Symbol
from sympy.matrices import Matrix, eye, zeros, ones, diag
from sympy import pretty_print

from SymbolicCollisions.core.cm_symbols import *
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator
from SymbolicCollisions.core.printers import *

matrixGenerator = MatrixGenerator(ex_D2Q9, ey_D2Q9, None, moments_dict[f'D2Q9'])
Mraw = matrixGenerator.get_raw_moments_matrix()

Nraw = matrixGenerator.get_shift_matrix()

Traw = matrixGenerator.get_raw_x_shift_moments_matrix()

Nraw_alternative = Traw * Mraw.inv()
Nraw_alternative_simplified = Matrix([round_and_simplify(Nraw_alternative[i, :]) for i in range(len(ex_D2Q9))])
# print_as_vector(f1, 'f1')

pretty_print(Mraw)
pretty_print(Nraw)
pretty_print(Nraw_alternative)

