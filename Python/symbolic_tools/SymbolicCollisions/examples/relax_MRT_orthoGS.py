
from sympy.matrices import eye
from SymbolicCollisions.core.cm_symbols import omega_v, S_relax_MRT_GS, M_ortho_GS
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_u2, print_as_vector, get_print_symbols_in_indx_notation

print("\n\n=== PRETTY CODE: relax relax_MRT_orthoGS ===\n\n")

DF_in_str = 'f_in'  # symbol defining DF
mom_DF_str = 'm'

# eq 10.30 from The Lattice Boltzmann Method: Principles and Practice
# T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
print("CudaDeviceFunction void relax_MRT_orthoGS("
      "real_t %s[9], "
      "real_t tau, "
      "\n{"
      % DF_in_str)


print("real_t %s = 1./tau;" % omega_v)
print("\nreal_t %s[9]; \n" % mom_DF_str)

DF = get_print_symbols_in_indx_notation(print_symbol=DF_in_str)
m_DF = get_print_symbols_in_indx_notation(print_symbol=mom_DF_str)
m = M_ortho_GS * DF

print("\n//orthogonal moments from density-probability functions")
print("//[m00, energy, energy^2, "
      "x momentum flux, x energy flux, "
      "y momentum flux, y energy flux, "
      "stress tensor (diagonal), stress tensor (off-diagonal)]")
print_as_vector(m, outprint_symbol=mom_DF_str)

print("\n//collision in orthogonal moments space")
print_as_vector(S_relax_MRT_GS * m_DF, outprint_symbol=mom_DF_str)

print("\n//back to density-probability functions")
DF = M_ortho_GS.inv() * m_DF
print_as_vector(DF, outprint_symbol=DF_in_str)

print("\n}\n")
