
from sympy.matrices import eye
from sympy.printing import print_ccode
from SymbolicCollisions.core.cm_symbols import omega_v, omega_b, Mraw_D2Q9, NrawD2Q9, S_relax_hydro_D2Q9, rho
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.hardcoded_results import hardcoded_cm_eq_incompressible_D2Q9, \
      hardcoded_F_cm_He_hydro_LB_incompressible_D2Q9, hardcoded_F_cm_Guo_hydro_LB_incompressible_D2Q9


print("\n\n=== PRETTY CODE: relax and collide ===\n\n")

q=9

pop_in_str = 'f_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations
F_cm_str = 'F_cm'

# eq : (eye(9)-S)*cm + S*cm_eq + (eye(9)-S/2.)*force_in_cm_space
print("CudaDeviceFunction void relax_and_collide_CM_hydro("
      f"real_t {pop_in_str}[{q}], "
      "real_t tau, "
      "vector_t u, "
      "vector_t Fhydro, "
      f"real_t {rho}) "
      "\n {")

print_u2()
print("real_t %s = 1./tau;" % omega_v)
# print("real_t bulk_visc = 1./6. ;")
# print("real_t %s = 1./(3*bulk_visc + 0.5);" % sb) 
print("real_t %s = omega_bulk;" % omega_b)  # s_b = 1./(3*bulk_visc + 0.5)
print("")

print_ccode(get_m00(q, pop_in_str), assign_to='real_t m00')

print("\nreal_t %s[9]; real_t %s[9]; real_t %s[9];\n" % (temp_pop_str, cm_eq_pop_str, F_cm_str))
print("for (int i = 0; i < 9; i++) {\n\t"
      "%s[i] = %s[i];}" % (temp_pop_str, pop_in_str))

populations = get_print_symbols_in_indx_notation(print_symbol=pop_in_str)
temp_populations = get_print_symbols_in_indx_notation(print_symbol=temp_pop_str)
cm_eq = get_print_symbols_in_indx_notation(print_symbol=cm_eq_pop_str)
F_cm = get_print_symbols_in_indx_notation(print_symbol=F_cm_str)
m = Mraw_D2Q9 * temp_populations

print("\n//raw moments from density-probability functions")
print("//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, outprint_symbol=pop_in_str)

print("\n//central moments from raw moments")
cm = NrawD2Q9 * populations
print_as_vector(cm, outprint_symbol=temp_pop_str)

print("\n//collision in central moments space")
print("//calculate equilibrium distributions in cm space")
# print_as_vector_re(get_cm_vector_from_discrete_def(get_pop_eq_hydro), cm_eq_pop_str)
print_as_vector(hardcoded_cm_eq_incompressible_D2Q9, cm_eq_pop_str)  # save time
print("//calculate forces in cm space")
# print_as_vector(get_cm_vector_from_discrete_def(get_force_Guo_second_order), F_cm_str)
# print_as_vector(get_cm_vector_from_continuous_def(get_continuous_force_He_MB), F_cm_str)
print_as_vector(hardcoded_F_cm_He_hydro_LB_incompressible_D2Q9, F_cm_str)  # save time

print("//collide eq: (eye(9)-S)*cm + S*cm_eq + (eye(9)-S/2.)*force_in_cm_space")
# cm_after_collision = (eye(9) - S_relax) * temp_populations + S_relax * cm_eq + (eye(9) - S_relax / 2) * F_cm
cm_after_collision = (eye(9) - S_relax_hydro_D2Q9) * temp_populations + S_relax_hydro_D2Q9 * cm_eq + (eye(
    9) - S_relax_hydro_D2Q9 / 2) * hardcoded_F_cm_He_hydro_LB_incompressible_D2Q9  # FOI - use `hardcoded` to skip zero terms
print_as_vector(cm_after_collision, outprint_symbol=pop_in_str)

print("\n//back to raw moments")
m = NrawD2Q9.inv() * populations
print_as_vector(m, outprint_symbol=temp_pop_str)

print("\n//back to density-probability functions")
populations = Mraw_D2Q9.inv() * temp_populations
print_as_vector(populations, outprint_symbol=pop_in_str)

print("\n}\n")
