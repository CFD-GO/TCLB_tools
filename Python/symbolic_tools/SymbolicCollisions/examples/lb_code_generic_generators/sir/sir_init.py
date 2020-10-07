from sympy.matrices import eye

from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v

from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00, get_m00
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation, get_vector_of_eq_central_moments
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht, print_sigma_sir
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
Init_S = Symbol('Init_S', positive=True)  # number
Init_I = Symbol('Init_I', positive=True)  # number
Init_R = Symbol('Init_R', positive=True)  # number

s_m_eq = get_vector_of_eq_central_moments(Init_S, Sigma2asSymbol)
i_m_eq = get_vector_of_eq_central_moments(Init_I, Sigma2asSymbol)
r_m_eq = get_vector_of_eq_central_moments(Init_R, Sigma2asSymbol)


# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

### GENERATE CODE ###
print(f"CudaDeviceFunction void set_eq_SIR(real_t rho) \n{{")
print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")

print(f"\treal_t {Init_S} = Init_S_Fraction/{rho};")
print(f"\treal_t {Init_I} = Init_R_Fraction/{rho};")
print(f"\treal_t {Init_R} = Init_R_Fraction/{rho};")

print_sigma_sir()

for m_eq, pop_in_str in zip([s_m_eq, i_m_eq, r_m_eq], ['s', 'i', 'r']):
    print(f"\t //--- processing {pop_in_str} ---")
    populations = get_print_symbols_in_m_notation(moments_order, pop_in_str)
    eq_pop_str = pop_in_str + "_eq_"
    eq_populations = get_print_symbols_in_m_notation(moments_order, eq_pop_str)

    print("\n\t//equilibrium")
    print_as_vector(m_eq, outprint_symbol="real_t " + eq_pop_str, output_order_of_moments=moments_order)

    print("\n\t//back to density-probability functions")
    print_as_vector(Mraw.inv() * eq_populations, outprint_symbol=pop_in_str, output_order_of_moments=rmoments_order)
    print("\n\n")

print("}\n")
