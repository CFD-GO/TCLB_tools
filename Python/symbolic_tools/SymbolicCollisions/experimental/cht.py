from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.DiscreteCMTransforms import DiscreteCMTransforms, get_mom_vector_from_discrete_def
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.cm_symbols import \
    F2D, dzeta2D, u2D, rho

from SymbolicCollisions.core.cm_symbols import rho, moments_dict
import time

lattice = 'D3Q27'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
# ccmt = ContinuousCMTransforms(dzeta2D, u2D, F2D, rho)
start = time.process_time()

print('\n\n// === discrete m === \n ')

from SymbolicCollisions.core.cm_symbols import e_D3Q7
print("moments: first order (linear) velocity expansion.")

dcmt = DiscreteCMTransforms(e_D3Q7, u3D, F3D, rho)
pop_eq = get_mom_vector_from_discrete_def(lambda i: dcmt.get_gamma_first_order_cht(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict['D3Q7'],
                                          serial_run=True)
print_as_vector(pop_eq, 'pop_eq_first_order', raw_output=True)

print('\n\n// === continous cm === \n ')

# to calculate particular moment
row = moments_dict['D2Q9'][0]
moment = ccmt.get_cm(row, ccmt.get_cht_DF)
print_as_vector(Matrix([moment]), 'particular_moment')



print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, \
    rho, cs2_thermal
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho, cs2=cs2_thermal)
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict['D2Q9'],
                                           serial_run=False)
print_as_vector(cm_eq, 'cm_eq')

print('---- CHT ----')
cm_cht_eq = get_mom_vector_from_continuous_def(ccmt.get_cht_DF,
                                               continuous_transformation=ccmt.get_cm,
                                               moments_order=moments_dict[lattice],
                                               serial_run=False)
print_as_vector(cm_cht_eq, 'cm_cht_eq', raw_output=False)

cm_cht_eq = get_mom_vector_from_continuous_def(ccmt.get_cht_DF,
                                               continuous_transformation=ccmt.get_m,
                                               moments_order=moments_dict[lattice],
                                               serial_run=False)
print_as_vector(cm_cht_eq, 'm_cht_eq', raw_output=False)

print(f'\n\n Done in {time.process_time() - start} [s].')
