from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix

from SymbolicCollisions.core.ContinuousCMTransforms import *
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho, \
    F2D, dzeta2D, u2D, rho
from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, M_ortho_GS
from SymbolicCollisions.core.cm_symbols import moments_dict, ex_D2Q9, Mraw_D2Q9, NrawD2Q9

import time

start = time.process_time()

"Corrections for diagonal third order velocity moments: feee" \
"eq 48, 49 from:"
"Coupling lattice Boltzmann model for simulation of thermal flows on stanard lattices" \
"by Q. Li, K. H. Luo, Y.L.He., Y.J. Gao, W.Q Tao, 2012"

phi_x = Symbol("phi_x")
phi_y = Symbol("phi_y")
C_D2Q9 = Matrix([-1 / 9 * phi_x,
                 -phi_x / 36 + phi_y / 4,
                 -phi_x / 36 - phi_y / 4,  # suprisingly the article applies asymmetric correction in x and y direction
                 -phi_x / 36 + phi_y / 4,
                 -phi_x / 36 - phi_y / 4,
                 phi_x / 18,
                 phi_x / 18,
                 phi_x / 18,
                 phi_x / 18, ])

print_as_vector(C_D2Q9.transpose() * Mraw_D2Q9, outprint_symbol="pop_in_str")

#
# print('\n//population_eq -> cm_eq - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)
# # cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)
# print_as_vector(cm_eq, 'cm_eq')


# import re
# re.findall(r'\d+', 'hello 42 I\'m a 32 string 30')
# ['42', '32', '30']
# # This would also match 42 from bla42bla. If you only want numbers delimited by word boundaries (space, period, comma), you can use \b :
#
# re.findall(r'\b\d+\b', 'he33llo 42 I\'m a 32 string 30')
# ['42', '32', '30']
# To end up with a list of numbers instead of a list of strings:
#
# >>> [int(s) for s in re.findall(r'\b\d+\b', 'he33llo 42 I\'m a 32 string 30')]
# [42, 32, 30]

