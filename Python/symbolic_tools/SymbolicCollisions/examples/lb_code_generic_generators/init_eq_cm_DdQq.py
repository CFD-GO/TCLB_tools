from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht, print_as_vector, get_print_symbols_in_m_notation, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix

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
Mraw = get_raw_moments_matrix(ex, ey, ez)
Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)

# from sympy import pprint
# pprint(Mraw)
# pprint(Nraw)

pop_in_str = 'h'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations

# GENERATE CODE
# print(f"CudaDeviceFunction void set_eq(real_t {pop_in_str}[{q}], real_t Xeq, vector_t u) \n{{")
print(f"CudaDeviceFunction void SetEquilibriumHeat(real_t H, real_t rho, vector_t u) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_sigma_cht()
print_u2(d)

# print(f"\treal_t {temp_pop_str}[{q}];\n")
# populations = get_print_symbols_in_indx_notation(q, pop_in_str)
# temp_populations = get_print_symbols_in_indx_notation(q, temp_pop_str)

populations = get_print_symbols_in_m_notation(moments_dict[f'D{d}Q{q}'], pop_in_str)
temp_populations = get_print_symbols_in_m_notation(moments_dict[f'D{d}Q{q}'], temp_pop_str)
for p in temp_populations:
    print(f"\treal_t {p};")

print("\n\t//equilibrium in central moments space")
print_as_vector(hardcoded_cm_eq, outprint_symbol=pop_in_str, moments_order=moments_dict[f'D{d}Q{q}'])

print("\n\t//back to raw moments")
# print_as_vector(Nraw.inv() * populations, outprint_symbol=temp_pop_str, moments_order=moments_dict[f'D{d}Q{q}'])
print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=temp_pop_str, moments_order=moments_dict[f'D{d}Q{q}'])  # shortcut
# print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=temp_pop_str, moments_order=moments_dict[f'D{d}Q{q}'], raw_output=True)  # shortcut

print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, outprint_symbol=pop_in_str, moments_order=moments_dict[f'D{d}Q{q}'])

print("\n}\n")
