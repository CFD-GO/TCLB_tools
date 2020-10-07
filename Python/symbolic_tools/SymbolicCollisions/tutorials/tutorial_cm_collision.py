from sympy.matrices import eye
from sympy.printing import print_ccode
from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00, moments_dict
from SymbolicCollisions.core.cm_symbols import dynamic_import
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.MatrixGenerator import MatrixGenerator

# inspired by:
# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space

# SETUP
d = 2
q = 9
model = 'hydro'  # choose from '['hydro', 'ade']

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None


def get_s_relax_switcher(choice):
    s_relax_switcher = {
        'hydro': ("SymbolicCollisions.core.cm_symbols", f"S_relax_hydro_D{d}Q{q}"),
        'ade': ("SymbolicCollisions.core.cm_symbols", f"S_relax_ADE_D{d}Q{q}"),
    }
    which_model = s_relax_switcher.get(choice, lambda: "Invalid argument")
    return dynamic_import(*which_model)


S_Relax = get_s_relax_switcher(model)

hardcoded_cm_eq = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_cm_eq_compressible_D{d}Q{q}")
hardcoded_F_cm = dynamic_import("SymbolicCollisions.core.hardcoded_results", f"hardcoded_F_cm_pf_D{d}Q{q}")
from SymbolicCollisions.core.cm_symbols import Force_str as F_str

# ARRANGE STUFF
matrixGenerator = MatrixGenerator(ex, ey, ez, moments_dict[f'D{d}Q{q}'])
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix()

# from sympy import pprint
# pprint(Mraw)  # see what you have done
# pprint(Nraw)

pop_in_str = 'x_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations


# GENERATE CODE
def make_header(choice):
    model_switcher = {
        'hydro': f"CudaDeviceFunction void relax_and_collide_hydro_with_F(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u, vector_t {F_str}) \n{{",
        'ade': f"CudaDeviceFunction void relax_and_collide_ADE(real_t {pop_in_str}[{q}], real_t {omega_ade}, vector_t u) \n{{",
    }
    result = model_switcher.get(choice, lambda: "Invalid argument")
    print(result)


make_header(model)

print("\t//=== THIS IS AUTOMATICALLY GENERATED CODE ===")
print_u2(d)
print_ccode(get_m00(q, pop_in_str), assign_to=f'\treal_t {m00}')


def make_variables(choice):
    model_switcher = {
        'hydro': f"\n\treal_t {temp_pop_str}[{q}];\n",
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
cm_eq = get_print_symbols_in_indx_notation(q, cm_eq_pop_str)
F_cm = get_print_symbols_in_indx_notation(q, F_str)
m = Mraw * temp_populations

print("\n\t//raw moments from density-probability functions")
# print("\t//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, outprint_symbol=pop_in_str)

print("\n\t//central moments from raw moments")
cm = Nraw * populations
print_as_vector(cm, outprint_symbol=temp_pop_str)

print("\n\t//collision in central moments space")
print("\t//collide")


def make_collision(choice):
    model_switcher = {
        # Relax 2nd moments for hydro, SOI
        'hydro': (eye(q) - S_Relax) * temp_populations
                 + S_Relax * hardcoded_cm_eq
                 + (eye(q) - S_Relax / 2) * hardcoded_F_cm,
        # Relax 1st moments for ADE, SOI without force
        'ade': None  # TODO: write the collision for advection-diffusion equation
    }
    # Get the function from switcher dictionary
    cm_after_collision = model_switcher.get(choice, lambda: "Invalid argument")
    print_as_vector(cm_after_collision, outprint_symbol=pop_in_str)

make_collision(model)

print("\n\t//back to raw moments")
# TODO: write 2 lines of code to
#  - reverse-transform central moments to raw moments
#  - print the result as c code;

print("\n\t//back to density-probability functions")
# TODO: write 2 lines of code to
#  - reverse-transform raw moments to density-probability functions
#  - print the result as c code;

print("\n}\n")
