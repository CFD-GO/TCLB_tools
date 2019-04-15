from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat

from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def

from sympy import Symbol
from SymbolicCollisions.core.cm_symbols import dzeta2D, e_D2Q9, u2D, F2D, rho, moments_dict, NrawD2Q9, Mraw_D2Q9, M_ortho_GS
from SymbolicCollisions.core.printers import print_as_vector

import time

start = time.process_time()

lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
ccmt = ContinuousCMTransforms(dzeta2D, u2D, F2D, rho)

print('\n\n// === discrete moments === \n ')
print('\n//moments from definition: k_mn = sum( (e_ix)^m (e_iy)^n * fun_i)')
print('\n\n// === BOUNDARY CONDITIONS === \n ')
print("discrete raw moments: velocity bc")
mom_bc = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_velocity_bc(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(mom_bc, 'drm_velocity_bc')

print("\n\n discrete raw moments: pressure bc")
mom_bc = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_pressure_bc(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])

print_as_vector(mom_bc, 'drm_pressure_bc')

print("\n\n continuous raw moments: concentration bc")

mom_bc = get_mom_vector_from_continuous_def(ccmt.get_bc_bb_concentration_cht,
                                          continuous_transformation=ccmt.get_m,
                                          moments_order=moments_dict[lattice])

print_as_vector(mom_bc, 'crm_concentration_cht_bc')

print("\n\n continuous raw moments: heat flux cht bc")
mom_bc = get_mom_vector_from_continuous_def(ccmt.get_force_He_MB,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])

print_as_vector(mom_bc, 'He forcing scheme')

mom_bc = get_mom_vector_from_continuous_def(ccmt.get_bc_bb_heat_flux_cht,
                                            continuous_transformation=ccmt.get_cm,
                                            moments_order=moments_dict[lattice])

print_as_vector(mom_bc, 'crm_heat_flux_cht_bc', raw_output=True)
