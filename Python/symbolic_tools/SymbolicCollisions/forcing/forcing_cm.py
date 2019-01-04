from SymbolicCollisions.core.sym_col_fun import *
from SymbolicCollisions.core.printers import print_as_vector, print_ccode
import numpy as np
import time
from sympy.utilities.iterables import flatten
from sympy import pretty_print
from sympy import *
start = time.process_time()

print('// === welcome to central moments space! === \n ')
print('// === discrete central moments ===\n ')


# print('\n//F_cm_Guo_extended from discrete def')
# F_cm_Guo_extended = get_mom_vector_from_discrete_def(get_discrete_force_Guo, discrete_transform=get_discrete_m)
# print_as_vector(F_cm_Guo_extended, 'F_cm', regex=True)
# #
# #
# print('\n//M*F_cm_Guo_extended')
# F_cm_Guo_extended = get_mom_vector_from_shift_Mat(get_discrete_force_Guo, Mat=NrawD2Q9 * Mraw_D2Q9)
# print_as_vector(F_cm_Guo_extended, 'F_cm', regex=True)
#
# print('\n\n// === continuous central moments === \n ')
#
#
# print('\n//Force -> Force_cm - from continuous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
# F_cm = get_mom_vector_from_continuous_def(get_continuous_force_Guo, continuous_transformation=get_continuous_m)
# print_as_vector(F_cm, 'F_cm', regex=True)

print('\n//population_eq -> cm_eq - from continuous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# # cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continous_transformation=get_continuous_cm)
cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_m)
print_as_vector(cm_eq, 'cm_eq', regex=True)

print('\n//Force -> Force_cm - from continuous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
# F_cm = get_mom_vector_from_continuous_def(get_continuous_force_He_first_order_MB, continous_transformation=get_continous_cm)
F_cm = get_mom_vector_from_continuous_def(get_continuous_force_He_hydro_DF, continuous_transformation=get_continuous_cm)
# F_cm = get_mom_vector_from_continuous_def(get_continuous_force_He_MB, continuous_transformation=get_continuous_m)
print_as_vector(F_cm, 'F_cm', regex=True)


print('\n\n Done in %s [s].'
      % str(time.process_time() - start))
