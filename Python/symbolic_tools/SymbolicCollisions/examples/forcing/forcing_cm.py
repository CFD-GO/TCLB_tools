from SymbolicCollisions.core.printers import print_as_vector
from SymbolicCollisions.core.ContinuousCMTransforms import \
    ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict, NrawD2Q9, Mraw_D2Q9

ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho)
import time

start = time.process_time()

lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

print("\n--- FORCES ---")
print('// === welcome to central moments space! === \n ')
print('// === discrete central moments ===\n ')

print('\n//F_cm_He_discrete')
F_cm_He = get_mom_vector_from_discrete_def(dcmt.get_force_He,
                                           discrete_transform=dcmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(F_cm_He, 'F_cm')

print('\n//N*M*F_He')
NMF_cm_He = get_mom_vector_from_shift_mat(dcmt.get_force_He, mat=NrawD2Q9 * Mraw_D2Q9)
print_as_vector(NMF_cm_He, 'F_cm')

print('\n//F_cm_Guo')
F_cm_Guo = get_mom_vector_from_discrete_def(dcmt.get_force_Guo,
                                            discrete_transform=dcmt.get_cm,
                                            moments_order=moments_dict[lattice])
print_as_vector(F_cm_Guo, 'F_cm')

print('\n//N*M*F_cm_Guo_second_order ')
NMF_cm_Guo = get_mom_vector_from_shift_mat(dcmt.get_force_Guo, mat=NrawD2Q9 * Mraw_D2Q9)
print_as_vector(NMF_cm_Guo, 'F_cm')

print('\n\n// === continuous central moments === \n ')

print('\n//Force -> Force_cm - from continuous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_Guo,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm')

print('\n//Force -> Force_cm - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_He_MB,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm')

print('\n//Force -> Force_cm - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_Guo,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm')

print('\n//Force -> Force_cm - from continous definition: \n'
      'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
      'where fun = forceM(rho,u,x,y) *(x-ux)^m *(y-uy)^n *(z-uz)^o ')
F_cm = get_mom_vector_from_continuous_def(ccmt.get_force_He_hydro_DF,
                                          continuous_transformation=ccmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm')

print(f'\n\n Done in {time.process_time() - start} [s].')
