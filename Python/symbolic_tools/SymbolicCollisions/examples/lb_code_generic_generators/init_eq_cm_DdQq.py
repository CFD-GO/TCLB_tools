from SymbolicCollisions.core.cm_symbols import rho, Enthalpy
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht, print_as_vector, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator

# SETUP
d = 3
q = 27

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")

if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None

e = dynamic_import("SymbolicCollisions.core.cm_symbols", f"e_D{d}Q{q}")

# hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}")
hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}")

# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex, ey, ez, moments_dict[f'D{d}Q{q}'])
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()


# from sympy import pprint
# pprint(Mraw)
# pprint(Nraw)

pop_in_str = 'h'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations

# GENERATE CODE
print(f"CudaDeviceFunction void set_eq(real_t {pop_in_str}[{q}], real_t {Enthalpy}, real_t {rho}, vector_t u) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_sigma_cht()
print_u2(d)

print(f"\treal_t {temp_pop_str}[{q}];\n")
populations = get_print_symbols_in_indx_notation(q, pop_in_str)
temp_populations = get_print_symbols_in_indx_notation(q, temp_pop_str)


print("\n\t//equilibrium in central moments space")
print_as_vector(hardcoded_cm_eq, outprint_symbol=pop_in_str)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * populations, outprint_symbol=temp_pop_str)

print("\n\t//back to density-probability functions")
eq_populations = Mraw.inv() * temp_populations
print_as_vector(eq_populations, outprint_symbol=pop_in_str)

print("\n}\n")
