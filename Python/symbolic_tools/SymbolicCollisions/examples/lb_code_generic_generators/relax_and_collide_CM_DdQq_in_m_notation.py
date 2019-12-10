from sympy.matrices import eye

from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00
from SymbolicCollisions.core.cm_symbols import Force_str as F_str
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator
from sympy.matrices import Matrix
import numpy as np
import pandas as pd

# inspired by:
# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space

model = 'cht'  # choose from '['hydro_compressible', 'hydro_incompressible', 'ade', 'ade_with_f', 'cht']
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


def get_s_relax_switcher(choice):
    s_relax_switcher = {
        'hydro_compressible': ("SymbolicCollisions.core.cm_symbols", f"S_relax_hydro_D{d}Q{q}"),
        'hydro_incompressible': ("SymbolicCollisions.core.cm_symbols", f"S_relax_hydro_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
        'cht': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
    }
    which_model = s_relax_switcher.get(choice, lambda: "Invalid argument")
    return dynamic_import(*which_model)


S_Relax = get_s_relax_switcher(model)

def get_cm_eq_and_F_cm_switcher(choice):
    cm_eq_switcher = {
        'hydro_compressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'hydro_incompressible': (
        "SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_incompressible_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'cht': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}"),
    }
    which_cm_eq = cm_eq_switcher.get(choice, lambda: "Invalid argument")
    hardcoded_cm_eq = dynamic_import(*which_cm_eq)

    F_cm_switcher = {
        'hydro_compressible': (
        "SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_hydro_compressible_D{d}Q{q}"),
        'hydro_incompressible': (
        "SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_hydro_incompressible_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_pf_D{d}Q{q}"),
        'ade': None,
        'cht': None,
    }
    which_F_cm = F_cm_switcher.get(choice, lambda: "Invalid argument")
    hardcoded_F_cm = \
        dynamic_import(*which_F_cm) if which_F_cm is not None \
            else Matrix(np.full((q, 1), 1, dtype=int))  # make dummy F

    # hardcoded_F_cm = dynamic_import(*which_F_cm)
    return hardcoded_cm_eq, hardcoded_F_cm


hardcoded_cm_eq, hardcoded_F_cm = get_cm_eq_and_F_cm_switcher(model)

# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

pop_in_str = 'g'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations


# GENERATE CODE
def make_header(choice):
    model_switcher = {
        'hydro_compressible': f"CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u, vector_t {F_str}) \n{{",
        'hydro_incompressible': f"CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u, vector_t {F_str}) \n{{",
        'ade_with_f': f"CudaDeviceFunction void relax_and_collide_ADE_with_F(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u, vector_t {F_str}) \n{{",
        'ade': f"CudaDeviceFunction void relax_and_collide_ADE(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u) \n{{",
        # 'cht': f"CudaDeviceFunction void relax_and_collide_ADE(real_t {pop_in_str}[<?R C( h_pop_size) ?>], real_t rho, real_t {omega_ade}, vector_t u) \n{{",
        'cht': f"CudaDeviceFunction void relax_and_collide_ADE_CM(real_t rho, real_t {omega_ade}, vector_t u) \n{{",
    }
    result = model_switcher.get(choice, lambda: "Invalid argument")
    print(result)


make_header(model)

print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
# print(f"real_t {sv} = omega;")
# print("real_t bulk_visc = 1./6. ;")
# print("real_t {sb} = 1./(3*bulk_visc + 0.5);")
# print(f"real_t {sb} = omega_bulk;\n")  # s_b = 1./(3*bulk_visc + 0.5)


if 'cht' in model:
    print_sigma_cht()

print_u2(d)

rtemp_populations = get_print_symbols_in_m_notation(rmoments_order, temp_pop_str)

populations = get_print_symbols_in_m_notation(moments_order, pop_in_str)
temp_populations = get_print_symbols_in_m_notation(moments_order, temp_pop_str)
cm_eq = get_print_symbols_in_m_notation(moments_order, cm_eq_pop_str)
F_cm = get_print_symbols_in_m_notation(moments_order, F_str)

if 'hydro_compressible' in model or 'hydro_incompressible' in model:
    print(f"\treal_t {omega_b} = 1.0;")

print(f"\treal_t {m00} = {sum(populations)};")

for t, p in zip(temp_populations, populations):
    print(f"\treal_t {t} = {p};")

print("\n\t//raw moments from density-probability functions")
print_as_vector(Mraw * rtemp_populations, outprint_symbol=pop_in_str, output_order_of_moments=moments_order)

print("\n\t//central moments from raw moments")
print_as_vector(Nraw * populations, outprint_symbol=temp_pop_str, output_order_of_moments=moments_order)

print("\n\t//collision in central moments space")
print("\t//collide")

def make_collision(choice):
    model_switcher = {
        # Relax 2nd moments for hydro, SOI
        'hydro_compressible': (eye(q) - S_Relax) * temp_populations
                              + S_Relax * hardcoded_cm_eq
                              + (eye(q) - S_Relax / 2) * hardcoded_F_cm,

        'hydro_incompressible': (eye(q) - S_Relax) * temp_populations
                                + S_Relax * hardcoded_cm_eq
                                + (eye(q) - S_Relax / 2) * hardcoded_F_cm,
        # Relax 1st moments for ADE, SOI
        'ade_with_f': (eye(q) - S_Relax) * temp_populations
                      + S_Relax * hardcoded_cm_eq
                      + (eye(q) - S_Relax / 2) * hardcoded_F_cm,
        # Relax 1st moments for ADE, SOI without force
        'ade': (eye(q) - S_Relax) * temp_populations
               + S_Relax * hardcoded_cm_eq,
        # Relax 1,3,5 moments for ADE, SOI without force
        'cht': (eye(q) - S_Relax) * temp_populations
               + S_Relax * hardcoded_cm_eq,
    }
    # Get the function from switcher dictionary
    cm_after_collision = model_switcher.get(choice, lambda: "Invalid argument")

    print_as_vector(cm_after_collision, outprint_symbol=pop_in_str, output_order_of_moments=moments_order)


make_collision(model)

print("\n\t//back to raw moments")
print_as_vector(Nraw.inv() * populations, outprint_symbol=temp_pop_str, output_order_of_moments=moments_order)

print("\n\t//back to density-probability functions")
print_as_vector(Mraw.inv() * temp_populations, outprint_symbol=pop_in_str, output_order_of_moments=rmoments_order)

print("\n}\n")
