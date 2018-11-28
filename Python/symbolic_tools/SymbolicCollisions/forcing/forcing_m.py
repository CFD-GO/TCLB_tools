from SymbolicCollisions.core.sym_col_fun import *
from SymbolicCollisions.core.printers import print_as_vector, print_ccode
from SymbolicCollisions.core.cm_symbols import Mraw, M_ortho_GS
import numpy as np
import time
from sympy.utilities.iterables import flatten
from sympy import pretty_print
from sympy import *
start = time.process_time()

print('// === welcome to moments space! === \n ')
print('// === discrete moments ===\n ')

T_raw_to_ortho = M_ortho_GS * Mraw.inv()

print('\n//F_m_Guo_extended')
F_m_Guo_extended = get_mom_vector_from_discrete_def(get_discrete_force_Guo_second_order, discrete_transform=get_discrete_m)
print_as_vector(F_m_Guo_extended, 'F_raw_m', regex=True)
print_as_vector(T_raw_to_ortho*F_m_Guo_extended, 'F_GS_m', regex=True)

print('\n//M*F_m_Guo_extended ')
F_m_Guo_extended = get_mom_vector_from_shift_Mat(get_discrete_force_Guo_second_order, Mat=Mraw)
print_as_vector(F_m_Guo_extended, 'F_raw_m', regex=True)
print_as_vector(T_raw_to_ortho*F_m_Guo_extended, 'F_GS_m', regex=True)


print('\n//M_ortho_GS*F_m_Guo_extended ')
F_m_Guo_extended = get_mom_vector_from_shift_Mat(get_discrete_force_Guo_second_order, Mat=M_ortho_GS)
print_as_vector(F_m_Guo_extended, 'F_GS_m', regex=True)


print('\n\n// === continuous moments === \n ')


print('\n//Gou`s forcing -> Force_m - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x)^m (y)^n ')
F_m = get_mom_vector_from_continuous_def(get_continuous_force_Guo, continuous_transformation=get_continuous_m)
print_as_vector(F_m, 'F_raw_m', regex=True)
print_as_vector(T_raw_to_ortho*F_m.transpose(), 'F_GS_m', regex=True)

print('\n//population_eq -> m_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fMB(rho,u,x,y) *(x)^m (y)^n ')
# cm_eq = get_cm_vector_from_continuous_def(get_continuous_Maxwellian_DF)
m_eq = get_mom_vector_from_continuous_def(get_continuous_hydro_DF, continuous_transformation=get_continuous_m)
print_as_vector(m_eq, 'm_raw_eq', regex=True)
print_as_vector(T_raw_to_ortho*m_eq.transpose(), 'm_GS_eq', regex=True)

print('\n//He`s forcing -> Force_m - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x)^m (y)^n ')
# F_cm = get_cm_vector_from_continuous_def(get_continuous_force_He_first_order_MB)
# F_cm = get_cm_vector_from_continuous_def(get_continuous_force_He_hydro_DF)
F_m = get_mom_vector_from_continuous_def(get_continuous_force_He_MB, continuous_transformation=get_continuous_m)
print_as_vector(F_m, 'F_raw_m', regex=True)
print_as_vector(T_raw_to_ortho*F_m.transpose(), 'F_GS_m', regex=True)

print('\n\n Done in %s [s].'
      % str(time.process_time() - start))
