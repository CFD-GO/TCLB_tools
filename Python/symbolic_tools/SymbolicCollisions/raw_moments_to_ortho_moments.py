
from SymbolicCollisions.core.sym_col_fun import *
from SymbolicCollisions.core.cm_symbols import *
from SymbolicCollisions.core.printers import print_as_vector

from sympy import pretty_print

print("\n\n=== from raw moments to ortho moments ===\n")
T_raw_to_ortho = M_ortho_GS * Mraw.inv()
# pretty_print(T_raw_to_ortho)

print("\n\n=== relax raw moments in ortho space ===\n")
S_relax2 = T_raw_to_ortho.inv() * S_relax_MRT_GS*T_raw_to_ortho
pretty_print(S_relax2)

print("\n\n deja-vu! \n")
pretty_print(S_relax)











print("\n\n=== PRETTY CODE: relax relax_MRT_relax_raw_mom_into_ortho ===\n\n")

DF_in_str = 'f_in'  # symbol defining DF
mom_DF_str = 'm'
mom_relaxed_DF_str = 'm_relaxed'

print("CudaDeviceFunction void relax_MRT_relax_raw_mom_into_ortho("
      "real_t %s[9], "
      "real_t tau, "
      "\n{"
      % DF_in_str)


print("\nreal_t %s = 1./tau;" % sv)
print("\nreal_t %s[9]; real_t %s[9]; \n" % (mom_DF_str, mom_relaxed_DF_str))

populations = get_DF(DF_in_str)
m_DF = get_DF(mom_DF_str)
m_relaxed_DF = get_DF(mom_relaxed_DF_str)
m = Mraw * populations

print("\n//raw moments from density-probability functions")
print("//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, print_symbol=mom_DF_str)


print("\n//collision in orthogonal moments space")
print_as_vector(S_relax2 * m_DF, print_symbol=mom_relaxed_DF_str)

print("\n//back to density-probability functions")
populations = Mraw.inv() * m_relaxed_DF
print_as_vector(populations, print_symbol=DF_in_str)

print("\n}\n")
