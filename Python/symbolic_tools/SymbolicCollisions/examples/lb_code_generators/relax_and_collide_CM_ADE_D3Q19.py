from sympy.matrices import eye

from SymbolicCollisions.core.cm_symbols import omega_v, omega_b, m00
from SymbolicCollisions.core.DiscreteCMTransforms import get_DF, get_m00

from SymbolicCollisions.core.printers import print_u2, print_as_vector, print_ccode

from SymbolicCollisions.core.MatrixGenerator import get_raw_moments_matrix, get_shift_matrix
from sympy import pprint

d = 3
q = 19

from SymbolicCollisions.core.cm_symbols import \
    ex_D3Q19 as ex, \
    ey_D3Q19 as ey, \
    ez_D3Q19 as ez

from SymbolicCollisions.core.cm_symbols import S_relax_ADE_D3Q19 as S_Relax_ADE

from SymbolicCollisions.core.hardcoded_results import \
    hardcoded_F_cm_pf_D3Q19 as hardcoded_F_cm

from SymbolicCollisions.core.hardcoded_results import \
    hardcoded_cm_eq_compressible_D3Q19 as hardcoded_cm_eq

Mraw = get_raw_moments_matrix(ex, ey, ez)
Nraw = get_shift_matrix(Mraw.inv(), ex, ey, ez)


forcing = False

# pprint(Mraw)
# pprint(Nraw)

print("\n\n=== PRETTY CODE: relax and collide ===\n\n")

pop_in_str = 'h_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations
F_cm_str = 'F_ade_cm'

# "Consistent Forcing Scheme in the cascaded LBM" L. Fei et al. 2017
# eqs 8-12 : (eye(q)-S)*cm + S*cm_eq + (eye(q)-S/2.)*force_in_cm_space


if forcing:
    print(f"CudaDeviceFunction void relax_and_collide_ADE_with_F(real_t {pop_in_str}[{q}], real_t omega, vector_t u, vector_t {F_cm_str}) \n{{")
else:
    print(f"CudaDeviceFunction void relax_and_collide_ADE(real_t {pop_in_str}[{q}], real_t {omega_v}, vector_t u) \n{{")


print_u2()
# print(f"real_t {sv} = omega;")
# print("real_t bulk_visc = 1./6. ;")
# print("real_t {sb} = 1./(3*bulk_visc + 0.5);")
# print(f"real_t {sb} = omega_bulk;\n")  # s_b = 1./(3*bulk_visc + 0.5)

print_ccode(get_m00(q, pop_in_str), assign_to=f'\treal_t {m00}')

if forcing:
    print(f"\n\treal_t {temp_pop_str}[{q}]; real_t {F_cm_str}[{q}];\n")
else:
    print(f"\n\treal_t {temp_pop_str}[{q}];\n")


print(f"\tfor (int i = 0; i < {q}; i++) {{\n\t"
      f"\t{temp_pop_str}[i] = {pop_in_str}[i];}}")

populations = get_DF(q, pop_in_str)
temp_populations = get_DF(q, temp_pop_str)
cm_eq = get_DF(q, cm_eq_pop_str)
F_cm = get_DF(q, F_cm_str)
m = Mraw * temp_populations

print("\n\t//raw moments from density-probability functions")
# print("\t//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, print_symbol=pop_in_str)
# print_as_vector(m, print_symbol=temp_pop_str, regex=False)

print("\n\t//central moments from raw moments")
cm = Nraw * populations
print_as_vector(cm, print_symbol=temp_pop_str)
# print_as_vector_old(cm, print_symbol=temp_pop_str, regex=False)

print("\n\t//collision in central moments space")
# print("//calculate equilibrium distributions in cm space")
# print("real_t {cm_eq_pop_str}[{q}];\n")
# print_as_vector(hardcoded_cm_eq, cm_eq_pop_str)  # save time, verbosity
# print("//calculate forces in cm space")
# print_as_vector(hardcoded_F_cm, F_cm_str)  # save time, verbosity
print("\t//collide")
if forcing:
    # Relax 1st moments for ADE, SOI
    cm_after_collision = (eye(q) - S_Relax_ADE) * temp_populations \
                         + S_Relax_ADE * hardcoded_cm_eq \
                         + (eye(q) - S_Relax_ADE / 2) * hardcoded_F_cm
else:
    cm_after_collision = (eye(q) - S_Relax_ADE) * temp_populations \
                         + S_Relax_ADE * hardcoded_cm_eq

print_as_vector(cm_after_collision, print_symbol=pop_in_str)

print("\n\t//back to raw moments")
m = Nraw.inv() * populations
print_as_vector(m, print_symbol=temp_pop_str)

print("\n\t//back to density-probability functions")
populations = Mraw.inv() * temp_populations
print_as_vector(populations, print_symbol=pop_in_str)

print("\n}\n")
