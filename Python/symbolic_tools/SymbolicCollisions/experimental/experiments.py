from SymbolicCollisions.core.sym_col_fun import *
from SymbolicCollisions.core.printers import print_as_vector, print_as_vector_new, print_ccode

import time
start = time.process_time()


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)
cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)

print_as_vector(cm_eq, 'cm_eq', regex=True)
print("\n-------------------------------------------------------")
print_as_vector_new(cm_eq, 'cm_eq', regex=True)
print("\n----------------------- without re -------------------------------")
print_as_vector_new(cm_eq, 'cm_eq', regex=False)


print('\n\n Done in %s [s].'
      % str(time.process_time() - start))


#
# print('\n//population_eq -> cm_eq - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)
# # cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)
# print_as_vector(cm_eq, 'cm_eq', regex=True)


# import re
# re.findall(r'\d+', 'hello 42 I\'m a 32 string 30')
# ['42', '32', '30']
# # This would also match 42 from bla42bla. If you only want numbers delimited by word boundaries (space, period, comma), you can use \b :
#
# re.findall(r'\b\d+\b', 'he33llo 42 I\'m a 32 string 30')
# ['42', '32', '30']
# To end up with a list of numbers instead of a list of strings:
#
# >>> [int(s) for s in re.findall(r'\b\d+\b', 'he33llo 42 I\'m a 32 string 30')]
# [42, 32, 30]