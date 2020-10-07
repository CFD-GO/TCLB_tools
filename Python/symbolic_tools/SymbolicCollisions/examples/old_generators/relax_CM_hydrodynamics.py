from sympy.printing import print_ccode
from SymbolicCollisions.core.cm_symbols import omega_v, omega_b, Mraw_D2Q9, NrawD2Q9, S_relax_hydro_D2Q9
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector, get_print_symbols_in_indx_notation
from SymbolicCollisions.core.hardcoded_results import hardcoded_cm_eq_incompressible_D2Q9, \
    hardcoded_F_cm_He_hydro_LB_incompressible_D2Q9

print("\n\n=== PRETTY CODE relax eq ===\n\n")
q = 9
pop_in_str = 'f_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations

# eq: -S*(cm - cm_eq)

print("CudaDeviceFunction void relax_CM_hydro("
      f"real_t {pop_in_str}[{q}], "
      "real_t tau, "
      # "vector_t Fhydro, "
      "vector_t u)"
      "\n{"
      )

print_u2()
print("real_t %s = 1./tau;" % omega_v)
# print("real_t bulk_visc = 1./6. ;")
# print("real_t %s = 1./(3*bulk_visc + 0.5);" % sb)  # s_b = 0.5; works good for some reason
print("real_t %s = omega_bulk;" % omega_b)  # s_b = 1./(3*bulk_visc + 0.5)
print("")

print_ccode(get_m00(q, pop_in_str), assign_to='real_t m00')

print("\nreal_t %s[9];" % temp_pop_str)
# print("\nreal_t %s[9];\n" % cm_eq_pop_str)

print("for (int i = 0; i < 9; i++) {\n\t"
      "%s[i] = %s[i];}" % (temp_pop_str, pop_in_str))

populations = get_print_symbols_in_indx_notation(print_symbol=pop_in_str)
temp_populations = get_print_symbols_in_indx_notation(print_symbol=temp_pop_str)
cm_eq = get_print_symbols_in_indx_notation(print_symbol=cm_eq_pop_str)
m = Mraw_D2Q9 * temp_populations

print("\n//raw moments from density-probability functions")
print("//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, outprint_symbol=pop_in_str)

print("\n//central moments from raw moments")
cm = NrawD2Q9 * populations
print_as_vector(cm, outprint_symbol=temp_pop_str)

print("\n//collision in central moments space")
print("//calculate equilibrium distributions in cm space")
# print_as_vector(get_cm_vector_from_discrete_def(get_pop_eq_hydro), cm_eq_pop_str)
print_as_vector(hardcoded_cm_eq_incompressible_D2Q9, cm_eq_pop_str)  # save time

print("//collide eq: -S*(cm - cm_eq)")
cm_after_collision = -S_relax_hydro_D2Q9 * (temp_populations - cm_eq)
# cm_after_collision = -S_relax * (temp_populations - hardcoded_cm_hydro_eq)
print_as_vector(cm_after_collision, outprint_symbol=pop_in_str)

print("\n//back to raw moments")
m = NrawD2Q9.inv() * populations
print_as_vector(m, outprint_symbol=temp_pop_str)

print("\n//back to density-probability functions")
populations = Mraw_D2Q9.inv() * temp_populations
print_as_vector(populations, outprint_symbol=pop_in_str)

print("\n}\n")
