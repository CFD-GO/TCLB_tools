from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho
from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, M_ortho_GS

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict
import time

from SymbolicCollisions.core.printers import get_print_symbols_in_indx_notation, get_print_symbols_in_m_notation
from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat


lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

start = time.process_time()

pop_in_str = 'f'
md = moments_dict['D2Q9']
mmd = Matrix(moments_dict['D2Q9'])

# populations_m = get_print_symbols_in_m_notation(e_D2Q9, pop_in_str)
populations_m = get_print_symbols_in_m_notation(mmd, pop_in_str)
populations_inx = get_print_symbols_in_indx_notation(e_D2Q9.shape[0], pop_in_str)
md = moments_dict['D2Q9']
mmd = Matrix(moments_dict['D2Q9'])
ed = [e_D2Q9[i, :] for i in range(e_D2Q9.shape[0])]

ed2 = get_print_symbols_in_m_notation(e_D2Q9, 'cantbempty')

print_as_vector(populations_m, outprint_symbol='g')
print()
print_as_vector(populations_m, outprint_symbol='g', output_order_of_moments=e_D2Q9)  # this is wrong!
print()
print_as_vector(populations_m, outprint_symbol='g', output_order_of_moments=mmd)
print('\n\n// ====== \n ')
print_as_vector(populations_inx)

print('\n\n// === continous cm === \n ')

print("\n--- EQUILIBRIA ---")

# to calculate particular moment
row = moments_dict['D2Q9'][0]
moment = ccmt.get_cm(row, ccmt.get_Maxwellian_DF)
print_as_vector(Matrix([moment]), 'particular_moment')


print('\n//population_eq -> cm_eq - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = fM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])

# cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
#                                            continuous_transformation=ccmt.get_cm,
#                                            moments_order=ed)

print_as_vector(cm_eq, 'cm_eq')
print("\n")
print_as_vector(cm_eq, 'cm_eq', output_order_of_moments=e_D2Q9)
print("\n")
print_as_vector(cm_eq, 'cm_eq', output_order_of_moments=mmd)

print(f'\n\n Done in {time.process_time() - start} [s].')

