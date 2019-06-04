from sympy.matrices import eye
from sympy.printing import print_ccode
from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00
from SymbolicCollisions.core.cm_symbols import Force_str as F_str
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF, get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix

# inspired by:
# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space

# SETUP
d = 3
q = 19

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None

Mraw = get_raw_moments_matrix(ex, ey, ez)
Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)

pop_in_str = 'x_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations

populations = get_DF(q, pop_in_str)
temp_populations = get_DF(q, temp_pop_str)

from sympy import pprint
pprint(Mraw)  # see what you have done
pprint(Nraw)
#
# x = Mraw.inv()
# pprint(x)  # see what you have done
# x = Nraw.inv()
# pprint(x)


print("\n\t//back to raw moments")
m = Nraw.inv() * populations
print_as_vector(m, print_symbol=temp_pop_str, raw_output=True)
#
print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, print_symbol=pop_in_str, raw_output=True)


print("\n\t//back to raw moments")
m = Nraw.inv() * populations
print_as_vector(m, print_symbol=temp_pop_str, raw_output=True)
#
print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, print_symbol=pop_in_str, raw_output=True)
