
import numpy as np
import pandas as pd
from SymbolicCollisions.core.cm_symbols import moments_dict, e_D3Q27, ex_D3Q27
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r

e_seed = [0, 1, -1]

grid = get_m_order_as_in_r(e_seed, e_seed, e_seed)
moments_order = grid.to_numpy()
print(f"order of moments: \n {moments_order}")
ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, e_D3Q27new = get_e_as_in_r(e_seed, e_seed, e_seed)
print(f"lattice velocities - e: \n {e_D3Q27new}")
print(e_D3Q27)

print(ex_D3Q27)
print(ex_D3Q27new)

print()
# > U
#    Var1 Var2 Var3
# 1     0    0    0
# 2     1    0    0
# 3    -1    0    0
# 4     0    1    0
# 5     1    1    0
# 6    -1    1    0
# 7     0   -1    0
# 8     1   -1    0
# 9    -1   -1    0
# 10    0    0    1
# 11    1    0    1
# 12   -1    0    1
# 13    0    1    1
# 14    1    1    1
# 15   -1    1    1
# 16    0   -1    1
# 17    1   -1    1
# 18   -1   -1    1
# 19    0    0   -1
# 20    1    0   -1
# 21   -1    0   -1
# 22    0    1   -1
# 23    1    1   -1
# 24   -1    1   -1
# 25    0   -1   -1
# 26    1   -1   -1
# 27   -1   -1   -1
