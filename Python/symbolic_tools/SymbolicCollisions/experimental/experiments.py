from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix

from SymbolicCollisions.core.ContinousCMTransforms import ContinousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho
from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, M_ortho_GS
from SymbolicCollisions.core.cm_symbols import moments_dict
import time

import time
start = time.process_time()

lattice = 'D2Q9'
ccmt = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])

# cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)

print_as_vector(cm_eq, 'cm_eq')
print("\n-------------------------------------------------------")
print_as_vector(cm_eq, 'cm_eq')
print("\n----------------------- without re -------------------------------")
print_as_vector(cm_eq, 'cm_eq')


print('\n\n Done in %s [s].'
      % str(time.process_time() - start))


#
# print('\n//population_eq -> cm_eq - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def(get_continuous_Maxwellian_DF, continuous_transformation=get_continuous_cm)
# # cm_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_cm)
# print_as_vector(cm_eq, 'cm_eq')


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