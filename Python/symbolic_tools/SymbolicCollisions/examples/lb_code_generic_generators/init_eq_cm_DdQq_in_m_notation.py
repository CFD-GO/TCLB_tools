from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht, print_as_vector, get_print_symbols_in_m_notation, get_e_as_in_r, get_m_order_as_in_r
from SymbolicCollisions.core.MatrixGenerator import get_e_as_in_r, get_m_order_as_in_r, MatrixGenerator
import numpy as np
import pandas as pd

# SETUP
m_seed = [0, 1, 2]
rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)
q, d = rmoments_order.shape
# hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results_m_notation", f"hardcoded_cm_eq_cht_D{d}Q{q}")

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

# from sympy import pprint
# pprint(Mraw)
# pprint(Nraw)

pop_in_str = 'h'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations

# GENERATE CODE
print(f"CudaDeviceFunction void SetEquilibriumHeat(real_t H, real_t rho, vector_t u) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_sigma_cht()
print_u2(d)

populations = get_print_symbols_in_m_notation(moments_dict[f'D{d}Q{q}'], pop_in_str)
temp_populations = get_print_symbols_in_m_notation(moments_dict[f'D{d}Q{q}'], temp_pop_str)

print("\n\t//equilibrium in central moments space")
hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}")
print_as_vector(hardcoded_cm_eq, outprint_symbol=pop_in_str, output_order_of_moments=moments_order)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * populations, outprint_symbol=f"real_t {temp_pop_str}", output_order_of_moments=moments_order)

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * temp_populations, outprint_symbol=pop_in_str, output_order_of_moments=rmoments_order)

print("\n}\n")
