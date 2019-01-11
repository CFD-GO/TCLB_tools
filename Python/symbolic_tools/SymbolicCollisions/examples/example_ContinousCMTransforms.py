from SymbolicCollisions.core.printers import print_as_vector

from SymbolicCollisions.core.hardcoded_results import \
    hardcoded_F_cm_Guo_hydro_LB_velocity_based_D2Q9, hardcoded_cm_eq_compressible_D2Q9, hardcoded_cm_eq_incompressible_D2Q9


from SymbolicCollisions.core.ContinousCMTransforms import ContinousCMTransforms, get_mom_vector_from_continuous_def_new
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, \
    F2D, dzeta2D, u2D, \
    rho, w, m00

from SymbolicCollisions.core.cm_symbols import moments_dict
import time


lattice = 'D2Q9'

# v1 = Matrix([3,2,1])
# v2 = Matrix([4,5,6])
#
# pprint(v1 - v2)
# pprint(v1.cross(v2))
# pprint(v1.dot(v2))


start = time.process_time()
# time.time() returns a value than can be interpreted as the current date and time.
# time.clock() returns a measure of how much CPU time has been used by the current process.


cm_i = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
# cm_i = ContinousCMTransforms(dzeta2D, u2D, F2D, rho)
# row = moments_dict['D2Q9'][0]
# test= cm_i.get_cm(row, cm_i.get_Maxwellian_DF)
# test=Matrix([test])
# print_as_vector(test, 'test', regex=True)

print("-----------------------------------------------------------------------------------------------------------")
# print('\n//population_eq -> cm_eq - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = fM(rho,u,x,y) *(x-ux)^m (y-uy)^n')
# cm_eq = get_mom_vector_from_continuous_def_new(cm_i.get_Maxwellian_DF,
#                                                continuous_transformation=cm_i.get_cm,
#                                                moments_order=moments_dict[lattice])


# # cm_eq = get_mom_vector_from_continuous_def_new(cm_i.get_hydro_DF,
# cm_eq = get_mom_vector_from_continuous_def_new(cm_i.get_force_He_MB,
# # cm_eq = get_mom_vector_from_continuous_def_new(cm_i.get_Maxwellian_DF,
#                                                continuous_transformation=cm_i.get_cm,
#                                                moments_order=moments_dict[lattice])
# print_as_vector(cm_eq, 'cm_eq', regex=True)

print("-----------------------------------------------------------------------------------------------------------")
# m_eq = get_mom_vector_from_continuous_def_new(cm_i.get_hydro_DF,
# m_eq = get_mom_vector_from_continuous_def_new(cm_i.get_Maxwellian_DF,
#                                                continuous_transformation=cm_i.get_m,
#                                                moments_order=moments_dict[lattice])
# print_as_vector(m_eq, 'm_eq', regex=True)

print("-----------------------------------------------------------------------------------------------------------")






# print_as_vector(cm_eq, 'cm_eq_no_regex', regex=False)

# print('\n//Force -> Force_cm - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
# F_cm = get_mom_vector_from_continuous_def_new(cm_i.get_force_He_MB,
F_cm = get_mom_vector_from_continuous_def_new(cm_i.get_force_He_hydro_DF,
                                              # F_cm=get_mom_vector_from_continuous_def_new(cm_i.get_force_Guo,
                                              continuous_transformation=cm_i.get_cm,
                                              moments_order=moments_dict[lattice])
print_as_vector(F_cm, 'F_cm', regex=True)






# print('\n//Force -> Force_cm - from continous definition: \n'
#       'k_mn = integrate(fun, (x, -oo, oo), (y, -oo, oo)) \n'
#       'where fun = forceM(rho,u,x,y) *(x-ux)^m (y-uy)^n ')
# F_cm_Guo_extended = get_mom_vector_from_discrete_def(get_discrete_force_Guo, discrete_transform=get_discrete_cm)
# print_as_vector(F_cm_Guo_extended, 'F_cm', regex=True)


# print('\n//N*M*F_He_continous ')
# from SymbolicCollisions.core.cm_symbols import \
#     ex_D3Q19 as ex, \
#     ey_D3Q19 as ey, \
#     ez_D3Q19 as ez

# from SymbolicCollisions.core.cm_symbols import \
#     ex_D2Q9 as ex, \
#     ey_D2Q9 as ey
# ez = None
# from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix
#
# Mraw = get_raw_moments_matrix(ex, ey, ez)
# Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)
#
# TO TYLKO DLA DYSKRETNEJ WERSJI!!!
# NMF_cm_He_original = get_mom_vector_from_shift_Mat(get_continuous_force_He_MB, Mat=Nraw * Mraw)
# print_as_vector(NMF_cm_He_original, 'F_cm', regex=True)  # produces looong expressions
#

print(f'\n\n Done in {time.process_time() - start} [s].')


