from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.cm_symbols import \
    F2D, dzeta2D, u2D, rho

from SymbolicCollisions.core.cm_symbols import rho, moments_dict
import time

lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
# ccmt = ContinuousCMTransforms(dzeta2D, u2D, F2D, rho)
start = time.process_time()

print('\n\n// === continous cm === \n ')

# to calculate particular moment
row = moments_dict['D2Q9'][0]
# moment = ccmt.get_cm(row, ccmt.get_Maxwellian_DF)
moment = ccmt.get_cm(row, ccmt.get_cht_DF)
print_as_vector(Matrix([moment]), 'particular_moment')



# print('\n//Force -> Force_cm - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
# F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_He_MB,
#                                           continuous_transformation=ccmt.get_cm,
#                                           moments_order=moments_dict[lattice])
# print_as_vector(F_cm, 'F_cm')


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

print(f'\n\n Done in {time.process_time() - start} [s].')
