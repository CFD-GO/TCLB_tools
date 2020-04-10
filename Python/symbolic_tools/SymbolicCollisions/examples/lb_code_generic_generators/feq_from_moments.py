from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import ux, uy, u2D, Fx, Fy, F2D, w_D2Q9
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator
from sympy.matrices import Matrix
import numpy as np
import pandas as pd
from sympy import Symbol

import time

start = time.process_time()

########################################

clip_z_dimension = True

m_seed = [0, 1, 2]
rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)

e_seed = [0, 1, -1]
ex, ey, ez, e_new = get_e_as_in_r(e_seed, e_seed, e_seed)

if clip_z_dimension:
    rmoments_order = rmoments_order[0:9]
    q, d = rmoments_order.shape
    d = 2
    ex = ex[0:9]
    ey = ey[0:9]
    ez = ez[0:9]
    e_D2Q9 = e_new[0:9, :]
else:
    q, d = rmoments_order.shape


print(f"order of rmoments: \n {pd.DataFrame.from_records(rmoments_order)}")
print(f"lattice velocities - e: \n {np.array(e_D2Q9)}")

######################
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)

print('\n\n// === continuous moments === \n ')
print('\n//population_eq -> m_eq - from continuous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fMB(rho,u,x,y) *(x)^m (y)^n ')

m_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                          continuous_transformation=ccmt.get_m,
                                          moments_order=rmoments_order)
print_as_vector(m_eq, 'm_raw_eq', raw_output=False, output_order_of_moments=rmoments_order)
print("--------------------------------------------------")

cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=rmoments_order)
print_as_vector(cm_eq, 'cm_eq', raw_output=False, output_order_of_moments=rmoments_order)

#########################
print("--------------------------------------------------")

matrixGenerator = MatrixGenerator(ex, ey, ez, rmoments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

feq = Mraw.inv() * m_eq.transpose()
print_as_vector(feq, outprint_symbol='f_eq_from_anal_mom', output_order_of_moments=rmoments_order)
print("--------------------------------------------------")

dcmt = DiscreteCMTransforms(e_D2Q9, Matrix([ux, uy, 0]), Matrix([Fx, Fy, 0]), rho, cs2=1./3., w=w_D2Q9)
discrete_edf = [dcmt.get_EDF(i) for i in range(0, 9)]
print_as_vector(Matrix([discrete_edf]), outprint_symbol=f"f_eq_2nd_order", output_order_of_moments=rmoments_order)
print("--------------------------------------------------")

f_eq_diff =feq -Matrix([discrete_edf]).transpose()
print_as_vector(f_eq_diff, outprint_symbol=f"feq_full_minus_f_eq_2nd_order", output_order_of_moments=rmoments_order)

print(f'\n\n Done in {time.process_time() - start} [s].')
print('bye')