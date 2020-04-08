from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import Symbol

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict
import time

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def


lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

start = time.process_time()
print('// === discrete (central) moments ===\n ')
print('// === welcome to TRT! === \n ')

print("moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'm_eq', output_order_of_moments=moments_dict[lattice])

print("TRT m antisymmetric - moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_TRT_antisymmetric(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'm_eq_antisymmetric', output_order_of_moments=moments_dict[lattice])

print("TRT m symmetric - moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_TRT_symmetric(i),
                                          discrete_transform=dcmt.get_m,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'm_eq_symmetric', output_order_of_moments=moments_dict[lattice])

print("=================================================================")

print("central moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma(i),
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'cm_eq', output_order_of_moments=moments_dict[lattice])

print("TRT cm antisymmetric - moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_TRT_antisymmetric(i),
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'cm_eq_antisymmetric', output_order_of_moments=moments_dict[lattice])

print("TRT cm symmetric - moments: second order (quadratic) velocity expansion.")
pop_eq = get_mom_vector_from_discrete_def(lambda i: Symbol('m00') * dcmt.get_gamma_TRT_symmetric(i),
                                          discrete_transform=dcmt.get_cm,
                                          moments_order=moments_dict[lattice])
print_as_vector(pop_eq, 'cm_eq_symmetric', output_order_of_moments=moments_dict[lattice])

print(f'\n\n Done in {time.process_time() - start} [s].')
