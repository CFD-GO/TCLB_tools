from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat

from sympy import Symbol
from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict, NrawD2Q9, Mraw_D2Q9
from SymbolicCollisions.core.printers import print_as_vector

import time

start = time.process_time()

lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)


print('\n\n// === continuous moments === \n ')
print('\n//population_eq -> m_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fMB(rho,u,x,y) *(x)^m (y)^n ')

m_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                          continuous_transformation=ccmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(m_eq, 'm_raw_eq')
T_raw_to_ortho = M_ortho_GS * Mraw_D2Q9.inv()
print_as_vector(T_raw_to_ortho*m_eq.transpose(), 'm_GS_eq')


print('// === welcome to cm! === \n ')


print(f'\n\n Done in {time.process_time() - start} [s].')
