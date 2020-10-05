from sympy.matrices import eye

from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00
from SymbolicCollisions.core.cm_symbols import Force_str as F_str
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator
from sympy.matrices import Matrix
import numpy as np
import pandas as pd


m_seed = [0, 1, 2]
rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)
q, d = rmoments_order.shape

moments_order = np.array(moments_dict[f'D{d}Q{q}'])
# print(f"order of moments | rmoments: \n "
#       f"{pd.concat([pd.DataFrame.from_records(moments_order),pd.DataFrame.from_records(rmoments_order)], axis=1)}")

e_seed = [0, 1, -1]
ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, e_D3Q27new = get_e_as_in_r(e_seed, e_seed, e_seed)
# print(f"lattice velocities - e: \n {np.array(e_D3Q27new)}")

matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
# Mraw = matrixGenerator.get_raw_moments_matrix()
# Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

pop_in_str = 'CS_'  # symbol defining populations
temp_pop_str = 'C_'  # symbol defining populations

populations = get_print_symbols_in_m_notation(moments_order, pop_in_str)
temp_populations = get_print_symbols_in_m_notation(moments_order, temp_pop_str)

#  const real_t CS_000 = C_000 ;
print(f"\treal_t H = {sum(populations)};")

for t, p in zip(temp_populations, populations):
    # print(f"\tconst real_t {p} = {t};")
    print(f"\tconst real_t {p} = 0;")
