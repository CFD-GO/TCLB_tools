from sympy.matrices import eye
from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00, moments_dict
from SymbolicCollisions.core.cm_symbols import Force_str as F_str
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator
from sympy.matrices import Matrix
import numpy as np

# inspired by:
# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space

# SETUP
d = 3
q = 27

pop_in_str = 'f'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations


model = 'hydro_incompressible'  # choose from '['hydro_compressible', 'hydro_incompressible', 'ade', 'ade_with_f', 'cht']

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None


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
        'hydro_incompressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_incompressible_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'cht': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_cht_D{d}Q{q}"),
    }
    which_cm_eq = cm_eq_switcher.get(choice, lambda: "Invalid argument")
    hardcoded_cm_eq = dynamic_import(*which_cm_eq)

    F_cm_switcher = {
        'hydro_compressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_hydro_compressible_D{d}Q{q}"),
        'hydro_incompressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_hydro_incompressible_D{d}Q{q}"),
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
matrixGenerator = MatrixGenerator(ex, ey, ez, moments_dict[f'D{d}Q{q}'])
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

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

print(f"\n\treal_t {temp_pop_str}[{q}];\n")

populations = get_print_symbols_in_indx_notation(q, pop_in_str)
temp_populations = get_print_symbols_in_indx_notation(q, temp_pop_str)
cm_eq = get_print_symbols_in_indx_notation(q, cm_eq_pop_str)
F_cm = get_print_symbols_in_indx_notation(q, F_str)

m = Mraw * temp_populations

print(f"\treal_t {m00} = {sum(populations)};")

print("\n\t//raw moments from density-probability functions")
print_as_vector(m, outprint_symbol=pop_in_str)

print("\n\t//central moments from raw moments")
cm = Nraw * populations
print_as_vector(cm, outprint_symbol=temp_pop_str)

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

    print_as_vector(cm_after_collision, outprint_symbol=pop_in_str)


print("\n\t//collision in central moments space")
print("\t//collide")
make_collision(model)

print("\n\t//back to raw moments")
m = Nraw.inv() * populations
print_as_vector(m, outprint_symbol=temp_pop_str)

print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, outprint_symbol=pop_in_str)

print("\n}\n")
