from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.cm_symbols import rho, moments_dict
import time


lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)

start = time.process_time()

print('\n\n// === continous cm === \n ')

# to calculate particular moment
row = moments_dict['D2Q9'][0]
moment = ccmt.get_cm(row, ccmt.get_Maxwellian_DF)
print_as_vector(Matrix([moment]), 'particular_moment')


print('\n//Force -> Force_cm - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_He_MB,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm')


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(cm_eq, 'cm_eq')


print(f'\n\n Done in {time.process_time() - start} [s].')
