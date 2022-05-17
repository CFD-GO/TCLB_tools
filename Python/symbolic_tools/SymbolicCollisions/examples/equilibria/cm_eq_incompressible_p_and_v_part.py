from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho
from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, M_ortho_GS

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict, S_relax_hydro_D2Q9
import time

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat


lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

start = time.process_time()

print('// === welcome to cm! === \n ')



print('\n//population_eq -> (central) m_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')

print('The raw moments for the pressure part')
m_eq_p = get_mom_vector_from_continuous_def(ccmt.get_incompressible_DF_part_p,
                                           continuous_transformation=ccmt.get_m,
                                           moments_order=moments_dict[lattice])

print_as_vector(m_eq_p, 'm_eq_p', output_order_of_moments=moments_dict[lattice])

print('The central moments for the velocity part')
cm_eq_gamma = get_mom_vector_from_continuous_def(ccmt.get_incompressible_DF_part_gamma,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])

print_as_vector(cm_eq_gamma, 'cm_eq_gamma', output_order_of_moments=moments_dict[lattice])

print('\n\n\n// let us relax the pressure part')
print('\n// pressure part - raw moments after relaxation')
print_as_vector(S_relax_hydro_D2Q9 * m_eq_p.transpose(), 'm_eq_p_after_collision', output_order_of_moments=moments_dict[lattice])


print('\n// pressure part - central moments')
cm_eq_p = get_mom_vector_from_continuous_def(ccmt.get_incompressible_DF_part_p,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(cm_eq_p, 'cm_eq_p', output_order_of_moments=moments_dict[lattice])
print('\n// pressure part - after relaxation after relaxation')
print_as_vector(S_relax_hydro_D2Q9 * cm_eq_p.transpose(), 'cm_eq_p_relaxation', output_order_of_moments=moments_dict[lattice])


print('\n\n\n// let us relax the velocity part')
print('\n// velocity part - central moments - after relaxation')
print_as_vector(S_relax_hydro_D2Q9 * cm_eq_gamma.transpose(), 'cm_eq_gamma_relaxed', output_order_of_moments=moments_dict[lattice])

print('\n\n\n// let us relax both parts together -  as we used to do in the article (eq C5 in appendix)')
print('\n// both parts - central moments - as we used to do in the article (eq C5 in appendix)')
cm_eq_both_parts = get_mom_vector_from_continuous_def(ccmt.get_incompressible_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])

print_as_vector(cm_eq_both_parts, 'cm_eq_both_parts', output_order_of_moments=moments_dict[lattice])
print('\n// cm_eq_both_parts - central moments after relaxation')
print_as_vector(S_relax_hydro_D2Q9 * cm_eq_both_parts.transpose(), 'cm_eq_both_parts_relaxation', output_order_of_moments=moments_dict[lattice])

print(f'\n\n Done in {time.process_time() - start} [s].')


print('\n\n\n// --------------- more EXPERIMENTS ------------------------')
print('\n// cm_eq_p_new')
cm_eq_p_new = get_mom_vector_from_discrete_def(lambda i: dcmt.get_EDF_p_with_u(i),
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(cm_eq_p_new, 'cm_eq_p_new', output_order_of_moments=moments_dict[lattice])

