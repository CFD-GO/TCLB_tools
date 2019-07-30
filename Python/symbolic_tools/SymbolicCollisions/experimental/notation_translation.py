from SymbolicCollisions.core.printers import print_as_vector
from sympy import Symbol
from sympy.matrices import Matrix
from SymbolicCollisions.core.printers import round_and_simplify
from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import rho
from SymbolicCollisions.core.cm_symbols import Temperature as T
from SymbolicCollisions.core.cm_symbols import cp, ux, uy

from SymbolicCollisions.core.cm_symbols import dzeta2D, e_D2Q9
from SymbolicCollisions.core.hardcoded_results import hardcoded_F_cm_hydro_density_based_D3Q19, \
    hardcoded_F_cm_Guo_hydro_LB_incompressible_D2Q9, \
    hardcoded_F_cm_hydro_density_based_D2Q9, \
    hardcoded_cm_eq_compressible_D2Q9,    hardcoded_cm_eq_compressible_D3Q19, \
    hardcoded_cm_eq_incompressible_D2Q9, \
    hardcoded_cm_eq_compressible_D2Q9_thermal


from SymbolicCollisions.core.cm_symbols import rho, moments_dict

from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF
from SymbolicCollisions.core.printers import print_u2, print_as_vector
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix
import re
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF

"""
Here are samples to facilitate debugging regex-printers ;)
"""

# SETUP
d = 2
q = 9
pop_in_str = 'x_in'  # symbol defining populations


# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None

# populations = get_DF(q, pop_in_str)
# # print_as_vector(populations, print_symbol='test')
#
# md = Matrix(moments_dict['D2Q9'])
#

md = e_D2Q9[4, :]

smd = [str(x) for x in md]
jsmd = ''.join(smd)
rjsmd = re.sub(r'-1', '2', jsmd)
stuff = Matrix([rjsmd])  # argh

print_as_vector(stuff, print_symbol='test')


# a =str(list(md._mat))
# pop_in_str = 'x_in'  # symbol defining populations
# temp_pop_str = 'temp'  # symbol defining populations
#
# populations = get_DF(q, pop_in_str)
# temp_populations = get_DF(q, temp_pop_str)
#
# Mraw = get_raw_moments_matrix(ex, ey, ez)
# m = Mraw * temp_populations
#
# # print("\n\t//raw moments from density-probability functions")
# # # print("\t//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
# print_as_vector(m, print_symbol=pop_in_str)

print("DONE")
