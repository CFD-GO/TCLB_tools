
from SymbolicCollisions.core.cm_symbols import Mraw_D2Q9, M_ortho_GS

from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict
from SymbolicCollisions.core.printers import print_as_vector
import time

start = time.process_time()

lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

T_raw_to_ortho = M_ortho_GS * Mraw_D2Q9.inv()

print('// === welcome to moments space! === \n ')
print('// === discrete moments ===\n ')

print('\n//F_m_Guo_extended')
F_m_Guo = get_mom_vector_from_discrete_def(dcmt.get_force_Guo,
                                           discrete_transform=dcmt.get_m,
                                           moments_order=moments_dict[lattice])
print_as_vector(F_m_Guo, 'F_cm')
print_as_vector(T_raw_to_ortho * F_m_Guo, 'F_GS_m')

print('\n//M*F_m_Guo_extended ')
MF_m_Guo = get_mom_vector_from_shift_mat(dcmt.get_force_Guo, mat=Mraw_D2Q9)
print_as_vector(MF_m_Guo, 'F_raw_m')
print_as_vector(T_raw_to_ortho * MF_m_Guo, 'F_GS_m')

print('\n//M_ortho_GS*F_m_Guo_extended ')
MF_m_Guo = get_mom_vector_from_shift_mat(dcmt.get_force_Guo, mat=M_ortho_GS)
print_as_vector(MF_m_Guo, 'F_GS_m')

print('\n\n// === continuous moments === \n ')

print('\n//Gou`s forcing -> Force_m - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x)^m (y)^n ')
F_m = get_mom_vector_from_continuous_def(ccmt.get_force_Guo,
                                         continuous_transformation=ccmt.get_m,
                                         moments_order=moments_dict[lattice])
print_as_vector(F_m, 'F_raw_m')
print_as_vector(T_raw_to_ortho * F_m.transpose(), 'F_GS_m')

print('\n//He`s forcing -> Force_m - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x)^m (y)^n ')
F_m = get_mom_vector_from_continuous_def(ccmt.get_force_He_MB,
                                         continuous_transformation=ccmt.get_m,
                                         moments_order=moments_dict[lattice])
print_as_vector(F_m, 'F_raw_m')
print_as_vector(T_raw_to_ortho * F_m.transpose(), 'F_GS_m')

print(f'\n\n Done in {time.process_time() - start} [s].')
