from sympy.matrices import eye
from sympy.printing import print_ccode
from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00, moments_dict
from SymbolicCollisions.core.cm_symbols import Force_str as F_str
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator
from sympy.matrices import Matrix
import numpy as np

# inspired by:
# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space

# SETUP
d = 3
q = 27
model = 'ade'  # choose from '['hydro_compressible', 'hydro_incompressible', 'ade', 'ade_with_f']

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None


def get_s_relax_switcher(choice):
    s_relax_switcher = {
        'hydro_compressible':   ("SymbolicCollisions.core.cm_symbols", f"S_relax_hydro_D{d}Q{q}"),
        'hydro_incompressible': ("SymbolicCollisions.core.cm_symbols", f"S_relax_hydro_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
    }
    which_model = s_relax_switcher.get(choice, lambda: "Invalid argument")
    return dynamic_import(*which_model)


S_Relax = get_s_relax_switcher(model)


def get_m_eq_and_F_m_switcher(choice):
    cm_eq_switcher = {
        'hydro_compressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'hydro_incompressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_incompressible_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_m_eq_D{d}Q{q}"),
    }
    which_m_eq = cm_eq_switcher.get(choice, lambda: "Invalid argument")
    hardcoded_m_eq = dynamic_import(*which_m_eq)


    F_m_switcher = {
        'hydro_compressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_hydro_density_based_D{d}Q{q}"),
        'hydro_incompressible': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_hydro_velocity_based_D{d}Q{q}"),
        'ade_with_f': ("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_pf_D{d}Q{q}"),
        'ade': None,
    }
    which_F_m = F_m_switcher.get(choice, lambda: "Invalid argument")
    hardcoded_F_m = \
        dynamic_import(*which_F_m) if which_F_m is not None \
        else Matrix(np.full((q, 1), 1, dtype=int))  # make dummy F


    # hardcoded_F_cm = dynamic_import(*which_F_cm)
    return hardcoded_m_eq, hardcoded_F_m


hardcoded_m_eq, hardcoded_F_m = get_m_eq_and_F_m_switcher(model)
# hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}")
# hardcoded_F_cm = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_pf_D{d}Q{q}")

# hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_incompressible_D{d}Q{q}")
# hardcoded_F_cm = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_He_hydro_LB_incompressible_D{d}Q{q}")


# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex, ey, ez, moments_dict[f'D{d}Q{q}'])
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

pop_in_str = 'x_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
m_eq_pop_str = 'm_eq'  # symbol defining populations


# GENERATE CODE
def make_header(choice):
    model_switcher = {
        'hydro_compressible': f"CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u, vector_t {F_str}) \n{{",
        'hydro_incompressible': f"CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u, vector_t {F_str}) \n{{",
        'ade_with_f': f"CudaDeviceFunction void relax_and_collide_ADE_with_F(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u, vector_t {F_str}) \n{{",
        'ade': f"CudaDeviceFunction void relax_and_collide_ADE(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u) \n{{",
    }
    result = model_switcher.get(choice, lambda: "Invalid argument")
    print(result)


make_header(model)

print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
# print(f"real_t {sv} = omega;")
# print("real_t bulk_visc = 1./6. ;")
# print("real_t {sb} = 1./(3*bulk_visc + 0.5);")
# print(f"real_t {sb} = omega_bulk;\n")  # s_b = 1./(3*bulk_visc + 0.5)
print_u2(d)
print_ccode(get_m00(q, pop_in_str), assign_to=f'\treal_t {m00}')


def make_variables(choice):
    model_switcher = {
        'hydro_compressible': f"\n\treal_t {temp_pop_str}[{q}];\n",
        'hydro_incompressible': f"\n\treal_t {temp_pop_str}[{q}];\n",
        'ade_with_f': f"\n\treal_t {temp_pop_str}[{q}];\n",
        'ade': f"\n\treal_t {temp_pop_str}[{q}];\n",
    }
    # Get the function from switcher dictionary
    result = model_switcher.get(choice, lambda: "Invalid argument")
    print(result)


make_variables(model)

print(f"\tfor (int i = 0; i < {q}; i++) {{\n\t"
      f"\t{temp_pop_str}[i] = {pop_in_str}[i];}}")

populations = get_print_symbols_in_indx_notation(q, pop_in_str)
temp_populations = get_print_symbols_in_indx_notation(q, temp_pop_str)
cm_eq = get_print_symbols_in_indx_notation(q, m_eq_pop_str)
F_cm = get_print_symbols_in_indx_notation(q, F_str)
m = Mraw * temp_populations

print("\n\t//raw moments from density-probability functions")
# print("\t//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, outprint_symbol=pop_in_str)


print("\n\t//collision in moments space")
# print("//calculate equilibrium distributions in cm space")
# print("real_t {cm_eq_pop_str}[{q}];\n")
# print_as_vector(hardcoded_cm_eq, cm_eq_pop_str)  # save time, verbosity
# print("//calculate forces in cm space")
# print("real_t {F_cm_str}[{q}];")
# print_as_vector(hardcoded_F_cm, F_cm_str)  # save time, verbosity
print("\t//collide")


def make_collision(choice):
    model_switcher = {
        # Relax 2nd moments for hydro, SOI
        'hydro_compressible': (eye(q) - S_Relax) * temp_populations
                              + S_Relax * hardcoded_m_eq
                              + (eye(q) - S_Relax / 2) * hardcoded_F_m,

        'hydro_incompressible': (eye(q) - S_Relax) * temp_populations
                                + S_Relax * hardcoded_m_eq
                                + (eye(q) - S_Relax / 2) * hardcoded_F_m,
        # Relax 1st moments for ADE, SOI
        'ade_with_f': (eye(q) - S_Relax) * temp_populations
                      + S_Relax * hardcoded_m_eq
                      + (eye(q) - S_Relax / 2) * hardcoded_F_m,
        # Relax 1st moments for ADE, SOI without force
        'ade': (eye(q) - S_Relax) * temp_populations
               + S_Relax * hardcoded_m_eq,
    }
    # Get the function from switcher dictionary
    cm_after_collision = model_switcher.get(choice, lambda: "Invalid argument")
    print_as_vector(cm_after_collision, outprint_symbol=pop_in_str)


make_collision(model)

print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, outprint_symbol=pop_in_str)

print("\n}\n")
