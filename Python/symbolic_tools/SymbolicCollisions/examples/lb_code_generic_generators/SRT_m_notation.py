from SymbolicCollisions.core.cm_symbols import m00
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict, omega_ade
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht, print_as_vector, get_print_symbols_in_m_notation
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator

import numpy as np

# SETUP
# d = 3
# q = 27

# # DYNAMIC IMPORTS
# ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
# ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
#
# if d == 3:
#     ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
# else:
#     ez = None
#
# e = dynamic_import("SymbolicCollisions.core.cm_symbols", f"e_D{d}Q{q}")

from SymbolicCollisions.core.printers import expand_grid_as_in_r, get_e_as_in_r

e_seed = [0, 1, -1]
m_seed = [0, 1, 2]

grid = expand_grid_as_in_r(m_seed, m_seed, m_seed)
moments_order = grid.to_numpy()
q, d = moments_order.shape
print(f"order of moments: \n {moments_order}")
# print(f"old order of moments: \n {np.array(moments_dict[f'D{d}Q{q}'])}")

ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, e_D3Q27new = get_e_as_in_r(e_seed, e_seed, e_seed)
print(f"lattice velocities - e: \n {np.array(e_D3Q27new)}")
# print(f"old lattice velocities - e: \n {np.array(e)}")
hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results_m_notation", f"hardcoded_cm_eq_cht_D{d}Q{q}")

# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix(Mraw.inv())

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
# for p in temp_populations:
#     print(f"\treal_t {p};")
#
# for peq in eq_populations:
#     print(f"\treal_t {peq};")

print("\n\t//equilibrium in central moments space")
# print_as_vector(hardcoded_cm_eq, outprint_symbol=pop_in_str, moments_order=moments_dict[f'D{d}Q{q}'])
print_as_vector(hardcoded_cm_eq, outprint_symbol=f"real_t {pop_eq_str}", moments_order=moments_order)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * eq_populations, outprint_symbol=f"real_t {temp_pop_str}", moments_order=moments_order)
# print_as_vector(Nraw.inv() * hardcoded_cm_eq, outprint_symbol=temp_pop_str, moments_order=moments_order)  # shortcut

print("\n\t//back to density-probability functions")
eq_populations = Mraw.inv() * temp_populations
print_as_vector(eq_populations, outprint_symbol=pop_eq_str, moments_order=moments_order)

print("\n\t//SRT collision")
for p, p_eq in zip(populations, eq_populations):
    print(f"\t{p} = {p} + {omega_ade}*({p_eq}-{p});")

print("\n}\n")
