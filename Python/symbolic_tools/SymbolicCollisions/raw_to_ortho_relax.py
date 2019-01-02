from SymbolicCollisions.core.cm_symbols import M_ortho_GS, Mraw_D2Q9, S_relax_MRT_GS
from sympy import pretty_print


print("\n\n=== is orthogonal and orthonormal? ===\n")
pretty_print(M_ortho_GS*M_ortho_GS.transpose())

print("\n\n=== from raw moments to ortho moments ===\n")
T_raw_to_ortho = M_ortho_GS * Mraw_D2Q9.inv()
pretty_print(T_raw_to_ortho)

print("\n\n=== relax raw moments in ortho space and go back to raw moments ===\n")
S_relax_ortho = T_raw_to_ortho.inv() * S_relax_MRT_GS * T_raw_to_ortho
pretty_print(S_relax_ortho)


