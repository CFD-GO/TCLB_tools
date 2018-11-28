from SymbolicCollisions.core.sym_col_fun import *
from SymbolicCollisions.core.printers import print_as_vector, print_ccode

import time
start = time.process_time()


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)
cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)
print_as_vector(cm_eq, 'cm_eq', regex=True)



print('\n\n Done in %s [s].'
      % str(time.process_time() - start))
