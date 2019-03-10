from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat

from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.ContinousCMTransforms import ContinousCMTransforms, get_mom_vector_from_continuous_def

from sympy import Symbol
from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict, NrawD2Q9, Mraw_D2Q9, M_ortho_GS
from SymbolicCollisions.core.printers import print_as_vector

import time

start = time.process_time()

lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
ccmt = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)

print('\n\n// === discrete moments === \n ')
print("moments: first order (linear) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_first_order(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'pop_eq_first_order')

print("moments: second order (quadratic) velocity expansion.")
print('\n//population_eq -> m_eq - by definition: k_mn = sum( (e_ix)^m (e_iy)^n * population_eq_i)')
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'pop_eq')

print('\n\n// === continuous moments === \n ')
print('\n//population_eq -> m_eq - from continuous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fMB(rho,u,x,y) *(x)^m (y)^n ')

m_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                          continuous_transformation=ccmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(m_eq, 'm_raw_eq')


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
m_eq = get_mom_vector_from_continuous_def(ccmt.get_hydro_DF,
                                           continuous_transformation=ccmt.get_m,
                                           moments_order=moments_dict[lattice])

print_as_vector(m_eq, 'm_raw_eq')

print("GS orthogonalization")
T_raw_to_ortho = M_ortho_GS * Mraw_D2Q9.inv()
print_as_vector(T_raw_to_ortho*m_eq.transpose(), 'm_GS_eq')


print(f'\n\n Done in {time.process_time() - start} [s].')
