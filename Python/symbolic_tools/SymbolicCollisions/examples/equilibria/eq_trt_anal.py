from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.DiscreteCMTransforms import DiscreteCMTransforms, get_mom_vector_from_discrete_def
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import ux, uy, u2D, Fx, Fy, F2D, w_D2Q9, moments_dict
from SymbolicCollisions.core.printers import print_as_vector, print_as_vector_latex, get_print_symbols_in_m_notation, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator, get_reverse_direction_idx, get_reverse_indices
from sympy.matrices import Matrix
import numpy as np
import pandas as pd
from sympy import Symbol

import time

start = time.process_time()

########################################

# clip_z_dimension = True
#
# m_seed = [0, 1, 2]
# rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)
#
# e_seed = [0, 1, -1]
# ex, ey, ez, e_new = get_e_as_in_r(e_seed, e_seed, e_seed)
#
# if clip_z_dimension:
#     rmoments_order = rmoments_order[0:9]
#     q, d = rmoments_order.shape
#     d = 2
#     ex = ex[0:9]
#     ey = ey[0:9]
#     ez = ez[0:9]
#     e_D2Q9 = e_new[0:9, :]
# else:
#     q, d = rmoments_order.shape

# order of moments given below is easier to read
from SymbolicCollisions.core.cm_symbols import ex_D2Q9 as ex, ey_D2Q9 as ey, ez_D2Q9 as ez
ez_D2Q9 = Matrix([0, 0, 0, 0, 0, 0, 0, 0, 0])
e_D2Q9 = ex.col_insert(1, ey)
e_D2Q9 = e_D2Q9.col_insert(2, ez)

rmoments_order = moments_dict['D2Q9']

rev_i = get_reverse_indices(e_D2Q9)
print(f"reverse direction indices: {rev_i}")
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
# print_as_vector(m_eq, 'm_raw_eq', raw_output=False, output_order_of_moments=rmoments_order)
print_as_vector_latex(m_eq, 'k^{H,eq}', output_order_of_moments=rmoments_order)

cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=rmoments_order)
# print_as_vector(cm_eq, 'cm_eq', raw_output=False, output_order_of_moments=rmoments_order)
print_as_vector_latex(cm_eq, '\\tilde{k}^{H,eq}', output_order_of_moments=rmoments_order)

print("--------------------------------------------------")

matrixGenerator = MatrixGenerator(ex, ey, ez, rmoments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

feq = Mraw.inv() * m_eq.transpose()
print_as_vector(feq, outprint_symbol='f_eq_from_anal_mom', output_order_of_moments=rmoments_order)
# print_as_vector_latex(feq, 'h^{eq}', output_order_of_moments=rmoments_order)

print("--------------------------------------------------")

dcmt = DiscreteCMTransforms(e_D2Q9, u3D, None, None)
# pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_TRT_antisymmetric(i),
#                                           discrete_transform=dcmt.get_m,
#                                           moments_order=rmoments_order)
# print_as_vector(pop_eq, 'dcm_eq_antisymmetric', output_order_of_moments=rmoments_order)
#

print("TRT cm antisymmetric - moments: full velocity expansion.")
# notice that the moments of non-eq DF is splited into sym and antisymmetric moments in the same way.
# feq = get_print_symbols_in_m_notation(rmoments_order, "f")

feq_symm = lambda i: (feq[i] + feq[rev_i[i]]) / 2
feq_antisymm = lambda i: (feq[i] - feq[rev_i[i]]) / 2

pop_eq_full = get_mom_vector_from_discrete_def(feq_symm,
                                          discrete_transform=dcmt.get_m,
                                          moments_order=rmoments_order)
# print_as_vector(pop_eq_full, 'm_eq_full_symmetric', output_order_of_moments=rmoments_order)
print_as_vector_latex(pop_eq_full, 'k^{H,seq}', output_order_of_moments=rmoments_order)

pop_eq_full = get_mom_vector_from_discrete_def(feq_antisymm,
                                          discrete_transform=dcmt.get_m,
                                          moments_order=rmoments_order)
# print_as_vector(pop_eq_full, 'm_eq_full_antisymmetric', output_order_of_moments=rmoments_order)
print_as_vector_latex(pop_eq_full, 'k^{H,aeq}', output_order_of_moments=rmoments_order)

pop_eq_full = get_mom_vector_from_discrete_def(feq_symm,
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=rmoments_order)
# print_as_vector(pop_eq_full, 'cm_eq_full_symmetric', output_order_of_moments=rmoments_order)
print_as_vector_latex(pop_eq_full, '\\tilde{k}^{H,seq}', output_order_of_moments=rmoments_order)

pop_eq_full = get_mom_vector_from_discrete_def(feq_antisymm,
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=rmoments_order)
# print_as_vector(pop_eq_full, 'cm_eq_full_antisymmetric', output_order_of_moments=rmoments_order)
print_as_vector_latex(pop_eq_full, '\\tilde{k}^{H,aeq}', output_order_of_moments=rmoments_order)
