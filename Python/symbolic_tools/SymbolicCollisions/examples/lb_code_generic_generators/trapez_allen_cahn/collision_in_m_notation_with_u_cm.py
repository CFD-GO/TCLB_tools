from sympy.matrices import eye

from SymbolicCollisions.core.cm_symbols import omega_ade, cs2
from SymbolicCollisions.core.printers import print_u2
from SymbolicCollisions.core.cm_symbols import moments_dict
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation, \
    get_vector_of_eq_moments, get_vector_of_eq_central_moments
from SymbolicCollisions.core.printers import print_sigma_sir
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator
from sympy.matrices import Matrix, diag
import numpy as np
import pandas as pd

from SymbolicCollisions.core.cm_symbols import Sigma2asSymbol, rho
from sympy import Symbol

# https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model
### PRELIMINARIES ###
clip_z_dimension = True

m_seed = [0, 1, 2]
rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)

e_seed = [0, 1, -1]
ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, e_D3Q27new = get_e_as_in_r(e_seed, e_seed, e_seed)

if clip_z_dimension:
    rmoments_order = rmoments_order[0:9]
    q, d = rmoments_order.shape
    d = 2
    ex_D3Q27new = ex_D3Q27new[0:9]
    ey_D3Q27new = ey_D3Q27new[0:9]
    ez_D3Q27new = ez_D3Q27new[0:9]
    e_D3Q27new = e_D3Q27new[0:9, :]
else:
    q, d = rmoments_order.shape

moments_order = np.array(moments_dict[f'D{d}Q{q}'])

print(f"order of moments | rmoments: \n "
      f"{pd.concat([pd.DataFrame.from_records(moments_order),pd.DataFrame.from_records(rmoments_order)], axis=1)}")

print(f"lattice velocities - e: \n {np.array(e_D3Q27new)}")

### PREPARE ENVIROMENT ###
# Relaxation_matrix = diag(1, omega_ade, omega_ade, 1, 1, 1, omega_ade, omega_ade, 1)
# Relaxation_matrix = diag(omega_ade, omega_ade, omega_ade, omega_ade, omega_ade, omega_ade, omega_ade, omega_ade, omega_ade)
omega_even = Symbol('omega_even', positive=True)
Relaxation_matrix = diag(omega_even, omega_ade, omega_ade, omega_even, omega_even, omega_even, omega_ade, omega_ade, omega_even)
Q = Symbol('Q', positive=True)
tilde_phi = Symbol('tilde_phi', positive=True)  # number


# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)



### GENERATE CODE ###
print(f"CudaDeviceFunction void relax_and_collide_TRT_CM_SOI(real_t rho, real_t {omega_ade}) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
f = Symbol('f', positive=True)  # fraction

print("\t// \"Optimal Stability of Advection-Diffusion Lattice Boltzmann Models \n"
      "\t// with Two Relaxation Times for Positive/Negative Equilibrium\" \n"
      "\t// by I. Ginzburg, D. d’Humières, A. Kuzmin, 2010\n"
      "\t// real_t omega_even = (-1.0 + 2.0/omega_ade)/(2.0*magic_parameter + 1.0);\n"
      "\t// real_t omega_even = 2. - omega_ade;\n"
      "\treal_t omega_even = 1.;\n"
      )


print_u2(d)
print(f"\treal_t {Sigma2asSymbol} = {cs2};")
print(f"\treal_t {Q} = getQ();")
rf_populations = get_print_symbols_in_m_notation(rmoments_order, f)

print(f"\treal_t {tilde_phi} = {sum(get_print_symbols_in_m_notation(moments_order, f))};")

mf_str = 'm_f_'
cmf_str = 'cm_f_'
cmf_eq_str = 'cm_f_eq_'
cmq_eq_str = 'cm_q_eq_'
collided_populations_str = "cm_star_"

mf = get_print_symbols_in_m_notation(moments_order, mf_str)
cmf = get_print_symbols_in_m_notation(moments_order, cmf_str)
cmf_eq = get_print_symbols_in_m_notation(moments_order, cmf_eq_str)
cmq_eq = get_print_symbols_in_m_notation(moments_order, cmq_eq_str)
cm_star = get_print_symbols_in_m_notation(moments_order, collided_populations_str)


print("\n\t//raw moments from density-probability functions")
print_as_vector(Mraw * rf_populations, outprint_symbol="real_t " + mf_str, output_order_of_moments=moments_order)

# same but in different order ...
# rmatrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, rmoments_order)
# rMraw = rmatrixGenerator.get_raw_moments_matrix()
# print_as_vector(rMraw * rf_populations, outprint_symbol="real_t " + 'rm', output_order_of_moments=rmoments_order)
print("\n\t//central moments from raw moments")
print_as_vector(Nraw * mf, outprint_symbol="real_t " + cmf_str, output_order_of_moments=moments_order)


print("\n\t//cm equilibrium moments ")

print_as_vector(get_vector_of_eq_central_moments(tilde_phi, Sigma2asSymbol), outprint_symbol="real_t " + cmf_eq_str, output_order_of_moments=moments_order)
print()
print_as_vector(get_vector_of_eq_central_moments(Q, Sigma2asSymbol), outprint_symbol="real_t " + cmq_eq_str, output_order_of_moments=moments_order)

print("\n\t//collide")
cm_after_collision = (eye(q) - Relaxation_matrix) * cmf + Relaxation_matrix * cmf_eq + cmq_eq
print_as_vector(cm_after_collision, outprint_symbol="real_t " + collided_populations_str,
                output_order_of_moments=moments_order)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * cm_star, outprint_symbol=mf_str, output_order_of_moments=moments_order)

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * mf, outprint_symbol=f, output_order_of_moments=rmoments_order)

print("}\n")

