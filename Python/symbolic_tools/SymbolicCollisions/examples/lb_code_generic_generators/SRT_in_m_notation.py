
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict, omega_ade
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht, print_as_vector, get_print_symbols_in_m_notation
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator
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

hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}")

# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)
# pprint(Nraw)

pop_in_str = 'h'  # symbol defining populations
pop_eq_str = 'heq'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations

# GENERATE CODE
print(f"CudaDeviceFunction void relax_and_collide_ADE_SRT_from_cm_eq(real_t rho, real_t {omega_ade},  vector_t u) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_sigma_cht()
print_u2(d)


populations = get_print_symbols_in_m_notation(moments_order, pop_in_str)
temp_populations = get_print_symbols_in_m_notation(moments_order, temp_pop_str)
eq_populations = get_print_symbols_in_m_notation(moments_order, pop_eq_str)

print(f"\treal_t H = {sum(populations)};")
print("\n\t//equilibrium in central moments space")
print_as_vector(hardcoded_cm_eq, outprint_symbol=f"real_t {pop_eq_str}", output_order_of_moments=moments_order)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * eq_populations, outprint_symbol=f"real_t {temp_pop_str}", output_order_of_moments=moments_order)
# print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=temp_pop_str, moments_order=moments_order)  # shortcut

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * temp_populations, outprint_symbol=pop_eq_str, output_order_of_moments=rmoments_order)

print("\n\t//SRT collision")
for p, p_eq in zip(populations, eq_populations):
    print(f"\t{p} = {p} + {omega_ade}*({p_eq}-{p});")

print("\n}\n")
