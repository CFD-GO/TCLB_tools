from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat

from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def

from sympy import Symbol
from sympy import Matrix
from SymbolicCollisions.core.cm_symbols import dzeta2D, e_D2Q9, u2D, F2D, rho, moments_dict, NrawD2Q9, Mraw_D2Q9, M_ortho_GS
from SymbolicCollisions.core.printers import print_as_vector

import time

start = time.process_time()

lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
ccmt = ContinuousCMTransforms(dzeta2D, u2D, F2D, rho)



# u2D * u2D

# e_D2Q9[2,:]
# Conservative phase-field lattice Boltzmann model for interface tracking equation
# PHYSICAL REVIEW E 91, 063309 (2015)
# Martin Geier, Abbas Fakhari and Taehun Lee


print('\n\n// === discrete moments === \n ')
print('\n//moments from definition: k_mn = sum( (e_ix)^m (e_iy)^n * fun_i)')
print('\n\n// === BOUNDARY CONDITIONS === \n ')
print("discrete raw moments: separation flux")
# mom_bc = get_mom_vector_from_discrete_def(lambda i: Symbol('H') * dcmt.get_heat_flux_bc(i),
#                                           discrete_transform=dcmt.get_cm,
#                                           moments_order=moments_dict[lattice])
# print_as_vector(mom_bc, 'dcm_heat_flux_cht_bc', raw_output=False)

mom_flux = get_mom_vector_from_discrete_def(dcmt.get_separation_flux,
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(mom_flux, 'mom_sep_flux', output_order_of_moments=moments_dict[lattice])
