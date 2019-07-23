from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF
from SymbolicCollisions.core.printers import print_u2, print_as_vector
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix

# SETUP
d = 3
q = 7

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None

hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}")
# hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}")

# ARRANGE STUFF
Mraw = get_raw_moments_matrix(ex, ey, ez)
Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)

# from sympy import pprint
# pprint(Mraw)
# pprint(Nraw)

pop_in_str = 'x_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations

# GENERATE CODE
print(f"CudaDeviceFunction void set_eq(real_t {pop_in_str}[{q}], real_t Xeq, vector_t u) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_u2(d)
print(f"\treal_t {m00} = Xeq;")
print(f"\treal_t {temp_pop_str}[{q}];\n")

populations = get_DF(q, pop_in_str)
temp_populations = get_DF(q, temp_pop_str)

print("\n\t//equilibrium in central moments space")
print_as_vector(hardcoded_cm_eq, print_symbol=pop_in_str)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * populations, print_symbol=temp_pop_str)
# print_as_vector(Nraw.inv() * hardcoded_cm_eq, print_symbol=temp_pop_str)  # shortcut

print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, print_symbol=pop_in_str)

print("\n}\n")
