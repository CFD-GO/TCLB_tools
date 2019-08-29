from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho
from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, M_ortho_GS

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict
import time

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat


lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

start = time.process_time()

print('// === welcome to cm! === \n ')
# print('// === discrete cm ===\n ')
#
#
# print('\n//population_eq -> cm_eq - by definition: k_mn = sum( (e_ix-ux)^m (e_iy-uy)^n * population_eq_i)')
# print("moments: first order (linear) velocity expansion.")
# pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_first_order(i),
#                                           discrete_transform=dcmt.get_cm,
#                                           moments_order=moments_dict[lattice])
# print_as_vector(pop_eq, 'pop_eq_first_order')
#
# print("moments: second order (quadratic) velocity expansion.")
# pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma(i),
#                                           discrete_transform=dcmt.get_cm,
#                                           moments_order=moments_dict[lattice])
# print_as_vector(pop_eq, 'pop_eq')
#
#
# print('\n//population -> cm - by definition: k_mn = sum( (e_ix-ux)^m (e_iy-uy)^n * population_i)')
# pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('%s[%d]' % ('pop', i)),
#                                           discrete_transform=dcmt.get_cm,
#                                           moments_order=moments_dict[lattice])
# print_as_vector(pop_eq, 'pop_cm')
#
# print('\n//velocity based hydrodynamic model: population_eq_pf -> cm_eq_pf - by definition: '
#       '\n//k_mn = sum( (e_ix-ux)^m (e_iy-uy)^n * population_eq_pf_i)')
# pop_eq = get_mom_vector_from_discrete_def(dcmt.get_EDF_incompressible,
#                                           discrete_transform=dcmt.get_cm,
#                                           moments_order=moments_dict[lattice])
# print_as_vector(pop_eq, 'cm_eq_pf')

print('\n\n// === continous cm === \n ')

print("\n--- EQUILIBRIA ---")

# to calculate particular moment
row = moments_dict['D2Q9'][0]
moment = ccmt.get_cm(row, ccmt.get_Maxwellian_DF)
print_as_vector(Matrix([moment]), 'particular_moment')


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])

print_as_vector(cm_eq, 'cm_eq')
print_as_vector(cm_eq, 'cm_eq', output_order_of_moments=moments_dict[lattice])

# print('\n//population_eq -> cm_eq - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
# cm_eq = get_mom_vector_from_continuous_def(ccmt.get_incompressible_DF,
#                                            continuous_transformation=ccmt.get_cm,
#                                            moments_order=moments_dict[lattice])
#
# print_as_vector(cm_eq, 'cm_eq')


print(f'\n\n Done in {time.process_time() - start} [s].')
