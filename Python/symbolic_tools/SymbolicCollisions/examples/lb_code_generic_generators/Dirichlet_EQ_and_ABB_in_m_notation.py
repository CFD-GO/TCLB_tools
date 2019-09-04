from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation, print_sigma_cht, print_u2
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.MatrixGenerator import get_e_as_in_r, get_m_order_as_in_r, MatrixGenerator
import numpy as np
import pandas as pd

# SETUP
m_seed = [0, 1, 2]
rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)
q, d = rmoments_order.shape

moments_order = np.array(moments_dict[f'D{d}Q{q}'])
print(f"order of moments | rmoments: \n "
      f"{pd.concat([pd.DataFrame.from_records(moments_order),pd.DataFrame.from_records(rmoments_order)], axis=1)}")

e_seed = [0, 1, -1]
ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, e_D3Q27new = get_e_as_in_r(e_seed, e_seed, e_seed)
print(f"lattice velocities - e: \n {np.array(e_D3Q27new)}")

# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

pop_in_str = 'h'  # symbol defining populations
eq_pop_str = 'heq'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations

populations = get_print_symbols_in_m_notation(moments_order, pop_in_str)
eq_populations = get_print_symbols_in_m_notation(moments_order, eq_pop_str)
temp_populations = get_print_symbols_in_m_notation(moments_order, temp_pop_str)

hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}")

print("\t//=== DIRICHLET - ABB ===")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print("\treal_t H = rho*cp*T")
print_sigma_cht()
print_u2(d)
print("\n\t//equilibrium in central moments space")
# print_as_vector(hardcoded_cm_eq, outprint_symbol=f"real_t {eq_pop_str}", output_order_of_moments=moments_order)

print("\n\t//back from cm_eq to raw moments")
print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=f"real_t {temp_pop_str}", output_order_of_moments=moments_order)
# print_as_vector(Nraw.inv() * eq_populations, outprint_symbol=f"real_t {temp_pop_str}", output_order_of_moments=moments_order)

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * temp_populations, outprint_symbol=f"real_t {eq_pop_str}", output_order_of_moments=rmoments_order)

print("\n\t//anti bounce back")
for p, p_eq in zip(populations, eq_populations):
    print(f"\t{p} = {-p} + 2 * {p_eq};")

print("\n}\n")


print("\t//=== DIRICHLET - EQ ===")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_sigma_cht()
print_u2(d)

print("\n\t//equilibrium in central moments space")
# print_as_vector(hardcoded_cm_eq, outprint_symbol=f"real_t {eq_pop_str}", output_order_of_moments=moments_order)
print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=f"real_t {temp_pop_str}", output_order_of_moments=moments_order)

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * temp_populations, outprint_symbol=pop_in_str, output_order_of_moments=rmoments_order)
print("\n}\n")


print("\t//=== IMPOSE HEAT FLUX ===")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_sigma_cht()
print_u2(d)

print("\n\t//equilibrium in central moments space")
# print_as_vector(hardcoded_cm_eq, outprint_symbol=f"real_t {eq_pop_str}", output_order_of_moments=moments_order)
print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=f"real_t {temp_pop_str}", output_order_of_moments=moments_order)

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * temp_populations, outprint_symbol=pop_in_str, output_order_of_moments=rmoments_order)
print("\n}\n")
